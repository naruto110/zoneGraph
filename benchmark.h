#ifndef BENCHMARK_H_
#define BENCHMARK_H_

#include "global.h"
#include "gapbs/pvector.h"
#include "gapbs/bitmap.h"
#include "gapbs/timer.h"
#include "gapbs/sliding_queue.h"
#include "gapbs/platform_atomics.h"
#include "gapbs/util.h"
// #include <condition_variable>
// #include <functional>
// #include <future>
// #include <chrono>

using namespace std;

typedef float ScoreT;
const float kDamp = 0.85;

static inline int32_t QueryVertexCount()
{
    return static_cast<int32_t>(global_vectex_vec.size());
}

struct QueryPrefetchedBlock {
    Vertex *vertex = nullptr;
    int32_t blk_idx = -1;
};

static inline vector<QueryPrefetchedBlock> prefetch_high_degree_query_blocks(int32_t block_limit)
{
    vector<QueryPrefetchedBlock> prefetched;
    if (block_limit <= 0) {
        return prefetched;
    }

    vector<Vertex *> candidates;
    candidates.reserve(global_vectex_vec.size());
    for (Vertex *v : global_vectex_vec) {
        if (v != nullptr && v->degree > CACHE_BLK_SIZE && !v->blk_list.empty()) {
            candidates.push_back(v);
        }
    }
    sort(candidates.begin(), candidates.end(), [](const Vertex *a, const Vertex *b) {
        if (a->degree != b->degree) {
            return a->degree > b->degree;
        }
        return a->id < b->id;
    });

    double start_t = get_current_time();
    for (Vertex *v : candidates) {
        unique_lock<std::shared_timed_mutex> lock(v->mutex);
        for (int32_t blkidx = 0;
             blkidx < static_cast<int32_t>(v->blk_list.size()) &&
             static_cast<int32_t>(prefetched.size()) < block_limit;
             blkidx++) {
            if (v->blk_list[blkidx].buf != nullptr) {
                continue;
            }
            v->blk_list[blkidx].buf = new int32_t[PM_BLK_SIZE];
            readBlockFromZoneInto(zns.dev_fd,
                                  v->blk_list[blkidx].buf,
                                  PM_BLK_SIZE,
                                  v->blk_list[blkidx].global_wp);
            prefetched.push_back({v, blkidx});
        }
        if (static_cast<int32_t>(prefetched.size()) >= block_limit) {
            break;
        }
    }
    double end_t = get_current_time();
    cout << "QUERY_PREFETCH_BLOCKS=" << prefetched.size() << endl;
    cout << "QUERY_PREFETCH_SEC=" << (end_t - start_t) << endl;
    return prefetched;
}

static inline void release_query_prefetched_blocks(const vector<QueryPrefetchedBlock> &prefetched)
{
    double start_t = get_current_time();
    int64_t released = 0;
    for (const QueryPrefetchedBlock &entry : prefetched) {
        if (entry.vertex == nullptr || entry.blk_idx < 0 ||
            entry.blk_idx >= static_cast<int32_t>(entry.vertex->blk_list.size())) {
            continue;
        }
        unique_lock<std::shared_timed_mutex> lock(entry.vertex->mutex);
        pm_blk &blk = entry.vertex->blk_list[entry.blk_idx];
        if (blk.buf != nullptr) {
            delete[] blk.buf;
            blk.buf = nullptr;
            released++;
        }
    }
    double end_t = get_current_time();
    cout << "QUERY_PREFETCH_RELEASED=" << released << endl;
    cout << "QUERY_PREFETCH_RELEASE_SEC=" << (end_t - start_t) << endl;
}

// Query-side adjacency iterator.
// Query prefetch may install clean read-only blocks into Vertex::blk_list[...].buf
// after maintenance has finished. The iterator never marks ZNS blocks invalid
// and never bumps global_cache_blk.
template <typename Fn>
static inline void for_each_neighbor_readonly(Vertex *v, Fn &&fn)
{
    std::shared_lock<std::shared_timed_mutex> lock(v->mutex);
    if (v->degree == 0) {
        return;
    }

    if (v->blk_list.empty()) {
        for (int32_t i = 0; i < v->degree; i++) {
            if (v->cache_buf[i] >= 0) {
                fn(v->cache_buf[i]);
            }
        }
        return;
    }

    int32_t tmp_blk_buf[PM_BLK_SIZE];
    for (int32_t blkidx = 0; blkidx < static_cast<int32_t>(v->blk_list.size()); blkidx++) {
        int32_t *blk_ptr = v->blk_list[blkidx].buf;
        if (blk_ptr == nullptr) {
            readBlockFromZoneInto(zns.dev_fd, tmp_blk_buf, PM_BLK_SIZE,
                                  v->blk_list[blkidx].global_wp);
            blk_ptr = tmp_blk_buf;
        }
        for (int32_t j = 0; j < PM_BLK_SIZE; j++) {
            if (blk_ptr[j] >= 0) {
                fn(blk_ptr[j]);
            }
        }
    }
    for (const auto &entry : v->pending_insert_by_block) {
        for (int32_t neighbor : entry.second) {
            if (neighbor >= 0) {
                fn(neighbor);
            }
        }
    }
    for (int32_t neighbor : v->pending_insert_unplaced) {
        if (neighbor >= 0) {
            fn(neighbor);
        }
    }
}

static inline bool triangle_block_contains_sorted(const int32_t *blk_ptr,
                                                  int32_t local_end,
                                                  int32_t target)
{
    int32_t lo = 0;
    int32_t hi = local_end - 1;
    while (lo <= hi) {
        int32_t mid = lo + (hi - lo) / 2;
        int32_t mid_idx = mid;
        while (mid_idx >= lo && blk_ptr[mid_idx] < 0) {
            mid_idx--;
        }
        if (mid_idx < lo) {
            lo = mid + 1;
            continue;
        }
        if (blk_ptr[mid_idx] == target) {
            return true;
        }
        if (blk_ptr[mid_idx] < target) {
            lo = mid + 1;
        } else {
            hi = mid_idx - 1;
        }
    }
    return false;
}

template <typename Fn>
static inline void triangle_for_each_neighbor_greater_than(Vertex *v,
                                                          int32_t lower_bound,
                                                          Fn &&fn)
{
    std::shared_lock<std::shared_timed_mutex> lock(v->mutex);
    if (v->degree == 0) {
        return;
    }

    if (v->blk_list.empty()) {
        for (int32_t i = 0; i < v->degree; i++) {
            int32_t neighbor = v->cache_buf[i];
            if (neighbor > lower_bound) {
                fn(neighbor);
            }
        }
        return;
    }

    int32_t tmp_blk_buf[PM_BLK_SIZE];
    for (int32_t blkidx = 0; blkidx < static_cast<int32_t>(v->blk_list.size()); blkidx++) {
        pm_blk &blk = v->blk_list[blkidx];
        if (blk.valid_count == 0 || blk.max_val <= lower_bound) {
            continue;
        }

        int32_t *blk_ptr = blk.buf;
        if (blk_ptr == nullptr) {
            readBlockFromZoneInto(zns.dev_fd, tmp_blk_buf, PM_BLK_SIZE, blk.global_wp);
            blk_ptr = tmp_blk_buf;
        }

        int32_t local_end = PM_BLK_SIZE;
        if (v->pma != nullptr) {
            local_end = pma_block_end(v->pma, blkidx) - pma_block_start(blkidx);
        }
        if (local_end <= 0) {
            continue;
        }

        if (blk.sorted) {
            int32_t lo = 0;
            int32_t hi = local_end;
            while (lo < hi) {
                int32_t mid = lo + (hi - lo) / 2;
                int32_t mid_idx = mid;
                while (mid_idx >= lo && blk_ptr[mid_idx] < 0) {
                    mid_idx--;
                }
                if (mid_idx < lo || blk_ptr[mid_idx] <= lower_bound) {
                    lo = mid + 1;
                } else {
                    hi = mid_idx;
                }
            }
            for (int32_t bit_idx = lo; bit_idx < local_end; bit_idx++) {
                int32_t neighbor = blk_ptr[bit_idx];
                if (neighbor < 0) {
                    continue;
                }
                if (neighbor <= lower_bound) {
                    continue;
                }
                fn(neighbor);
            }
            continue;
        }

        for (int32_t bit_idx = 0; bit_idx < local_end; bit_idx++) {
            int32_t neighbor = blk_ptr[bit_idx];
            if (neighbor > lower_bound) {
                fn(neighbor);
            }
        }
    }

    for (const auto &entry : v->pending_insert_by_block) {
        for (int32_t neighbor : entry.second) {
            if (neighbor > lower_bound) {
                fn(neighbor);
            }
        }
    }
    for (int32_t neighbor : v->pending_insert_unplaced) {
        if (neighbor > lower_bound) {
            fn(neighbor);
        }
    }
}

static inline bool triangle_contains_edge(Vertex *v, int32_t target)
{
    std::shared_lock<std::shared_timed_mutex> lock(v->mutex);
    if (v->degree == 0) {
        return false;
    }

    if (v->blk_list.empty()) {
        for (int32_t i = 0; i < v->degree; i++) {
            if (v->cache_buf[i] == target) {
                return true;
            }
        }
        return false;
    }

    int32_t tmp_blk_buf[PM_BLK_SIZE];
    for (int32_t blkidx = 0; blkidx < static_cast<int32_t>(v->blk_list.size()); blkidx++) {
        pm_blk &blk = v->blk_list[blkidx];
        if (blk.valid_count == 0 || target < blk.min_val || target > blk.max_val) {
            continue;
        }

        int32_t *blk_ptr = blk.buf;
        if (blk_ptr == nullptr) {
            readBlockFromZoneInto(zns.dev_fd, tmp_blk_buf, PM_BLK_SIZE, blk.global_wp);
            blk_ptr = tmp_blk_buf;
        }

        int32_t local_end = PM_BLK_SIZE;
        if (v->pma != nullptr) {
            local_end = pma_block_end(v->pma, blkidx) - pma_block_start(blkidx);
        }
        if (local_end <= 0) {
            continue;
        }

        if (blk.sorted) {
            if (triangle_block_contains_sorted(blk_ptr, local_end, target)) {
                return true;
            }
            continue;
        }

        for (int32_t bit_idx = 0; bit_idx < local_end; bit_idx++) {
            if (blk_ptr[bit_idx] == target) {
                return true;
            }
        }
    }

    for (const auto &entry : v->pending_insert_by_block) {
        for (int32_t neighbor : entry.second) {
            if (neighbor == target) {
                return true;
            }
        }
    }
    for (int32_t neighbor : v->pending_insert_unplaced) {
        if (neighbor == target) {
            return true;
        }
    }
    return false;
}

int64_t CountTrianglesById(void)
{
    unique_lock<std::shared_timed_mutex> graph_lock(global_query_update_mutex);
    const int32_t vertex_count = QueryVertexCount();
    int64_t total_triangles = 0;

    #pragma omp parallel for reduction(+ : total_triangles) schedule(dynamic, 64)
    for (int32_t u = 0; u < vertex_count; u++) {
        Vertex *u_vertex = global_vectex_vec[u];
        if (u_vertex == nullptr || u_vertex->degree < 2) {
            continue;
        }

        vector<int32_t> forward_neighbors;
        forward_neighbors.reserve(static_cast<size_t>(u_vertex->degree));
        triangle_for_each_neighbor_greater_than(u_vertex, u, [&](int32_t neighbor) {
            if (neighbor >= 0 && neighbor < vertex_count) {
                forward_neighbors.push_back(neighbor);
            }
        });

        if (forward_neighbors.size() < 2) {
            continue;
        }

        sort(forward_neighbors.begin(), forward_neighbors.end(),
             [](int32_t lhs, int32_t rhs) { return lhs < rhs; });
        forward_neighbors.erase(unique(forward_neighbors.begin(),
                                       forward_neighbors.end()),
                                forward_neighbors.end());

        for (size_t i = 0; i + 1 < forward_neighbors.size(); i++) {
            int32_t v = forward_neighbors[i];
            Vertex *v_vertex = global_vectex_vec[v];
            for (size_t j = i + 1; j < forward_neighbors.size(); j++) {
                int32_t w = forward_neighbors[j];
                if (triangle_contains_edge(v_vertex, w)) {
                    total_triangles++;
                }
            }
        }
    }

    return total_triangles;
}

// 在当前图存储之上直接执行 PageRank。
// 这里既支持从 cache_buf 读取，也支持从 zone 中按块读取邻接。
ScoreT * PageRankPullGS(int max_iters, double epsilon = 1e-4) 
{
    unique_lock<std::shared_timed_mutex> graph_lock(global_query_update_mutex);
    omp_set_num_threads(THREAD_NUM);
    const int32_t vertex_count = QueryVertexCount();
    const ScoreT init_score = 1.0f / vertex_count;
    const ScoreT base_score = (1.0f - kDamp) / vertex_count;
    ScoreT *scores;
    ScoreT *outgoing_contrib;
    scores = (ScoreT *) malloc(sizeof(ScoreT) * vertex_count);
    outgoing_contrib = (ScoreT *) malloc(sizeof(ScoreT) * vertex_count);

    #pragma omp parallel for
    // #pragma omp parallel for num_threads(THREAD_NUM)
    for (int32_t vidx = 0; vidx < vertex_count; vidx++) {
        scores[vidx] = init_score;
        outgoing_contrib[vidx] = init_score;
    }


    for (int iter=0; iter < max_iters; iter++) {
        double error = 0;
        #pragma omp parallel for reduction(+ : error) schedule(dynamic, THREAD_NUM)
        for (int32_t vidx=0; vidx < vertex_count; vidx++) {
        // printf("I am thread %d / %d \n", omp_get_thread_num(),omp_get_num_threads());
        // cout << "thread: " << omp_get_thread_num() << " total threads: " << omp_get_num_threads();
            ScoreT incoming_total = 0;
            Vertex *v = global_vectex_vec[vidx];
            int32_t degree = 0;
            for_each_neighbor_readonly(v, [&](int32_t neighbor) {
                incoming_total += outgoing_contrib[neighbor];
                degree++;
            });

            ScoreT old_score = scores[vidx];
            scores[vidx] = base_score + kDamp * incoming_total;
            error += fabs(scores[vidx] - old_score);
            outgoing_contrib[vidx] = degree == 0 ? 0 : scores[vidx] / degree;
            // delete[] tmp_blk_ptr;
        }
        // cout << "iter: " << iter << " total_num_access: " << total_num_access << endl;
    //  printf(" %2d    %lf\n", iter, error);
    //  if (error < epsilon)
    //    break;
    }

    free(outgoing_contrib);
    return scores;
}



// Bottom-up BFS 单步。
// front 是当前前沿，next 保存下一轮被激活的点。
int32_t BUStep(pvector<int32_t> &parent, Bitmap &front, Bitmap &next) {
    int32_t awake_count = 0;
    const int32_t vertex_count = QueryVertexCount();
    next.reset();
    #pragma omp parallel for reduction(+ : awake_count) schedule(dynamic, 1024)
    for (int32_t uidx=0; uidx < vertex_count; uidx++) {
        if (parent[uidx] < 0) {
            Vertex *v = global_vectex_vec[uidx];
            bool found = false;
            for_each_neighbor_readonly(v, [&](int32_t neighbor) {
                if (!found && front.get_bit(neighbor)) {
                    parent[uidx] = neighbor;
                    awake_count++;
                    next.set_bit(uidx);
                    found = true;
                }
            });
        }
    }
  return awake_count;
}


int64_t TDStep(pvector<int32_t> &parent, SlidingQueue<int32_t> &queue) {
    int64_t scout_count = 0;
    #pragma omp parallel
    {
        QueueBuffer<int32_t> lqueue(queue);
        #pragma omp for reduction(+ : scout_count)
        for (auto q_iter = queue.begin(); q_iter < queue.end(); q_iter++) {
            int32_t uidx = *q_iter;
            Vertex *v = global_vectex_vec[uidx];
            for_each_neighbor_readonly(v, [&](int32_t neighbor) {
                    int32_t curr_val = parent[neighbor];
                    if (curr_val < 0) {
                        if (compare_and_swap(parent[neighbor], curr_val, uidx)) {
                            lqueue.push_back(neighbor);
                            scout_count += -curr_val;
                        }
                    }
            });
        }
        lqueue.flush();
    }
    return scout_count;
}


void QueueToBitmap(const SlidingQueue<int32_t> &queue, Bitmap &bm) {
    #pragma omp parallel for
    for (auto q_iter = queue.begin(); q_iter < queue.end(); q_iter++) {
        int32_t u = *q_iter;
        bm.set_bit_atomic(u);
    }
}

void BitmapToQueue(const Bitmap &bm, SlidingQueue<int32_t> &queue) {
    const int32_t vertex_count = QueryVertexCount();
    #pragma omp parallel
    {
        QueueBuffer<int32_t> lqueue(queue);
        #pragma omp for
        for (int32_t n=0; n < vertex_count; n++)
        if (bm.get_bit(n))
            lqueue.push_back(n);
        lqueue.flush();
    }
    queue.slide_window();
}

pvector<int32_t> InitParent(void) {
    const int32_t vertex_count = QueryVertexCount();
    pvector<int32_t> parent(vertex_count);
    #pragma omp parallel for
    for (int32_t n=0; n < vertex_count; n++) {
        std::shared_lock<std::shared_timed_mutex> lock(global_vectex_vec[n]->mutex);
        parent[n] = global_vectex_vec[n]->degree != 0 ? -global_vectex_vec[n]->degree : -1;
    }
    return parent;
}

pvector<int32_t> DOBFS(int32_t source, int alpha = 15, int beta = 18) {
//  PrintStep("Source", static_cast<int32_t>(source));
    unique_lock<std::shared_timed_mutex> graph_lock(global_query_update_mutex);
    const int32_t vertex_count = QueryVertexCount();
    Timer t;
    t.Start();
    pvector<int32_t> parent = InitParent();
    t.Stop();
    //  PrintStep("i", t.Seconds());
    parent[source] = source;
    SlidingQueue<int32_t> queue(vertex_count);
    queue.push_back(source);
    queue.slide_window();
    Bitmap curr(vertex_count);
    curr.reset();
    Bitmap front(vertex_count);
    front.reset();
    int64_t edges_to_check = global_edge_num;
    int64_t scout_count = 0;
    {
        std::shared_lock<std::shared_timed_mutex> lock(global_vectex_vec[source]->mutex);
        scout_count = global_vectex_vec[source]->degree;
    }
    while (!queue.empty()) {
        if (scout_count > edges_to_check / alpha) {
        int32_t awake_count, old_awake_count;
        TIME_OP(t, QueueToBitmap(queue, front));
    //      PrintStep("e", t.Seconds());
        awake_count = queue.size();
        queue.slide_window();
        do {
            t.Start();
            old_awake_count = awake_count;
            awake_count = BUStep(parent, front, curr);
            front.swap(curr);
            t.Stop();
    //        PrintStep("bu", t.Seconds(), awake_count);
        } while ((awake_count >= old_awake_count) ||
                (awake_count > vertex_count / beta));
        TIME_OP(t, BitmapToQueue(front, queue));
    //      PrintStep("c", t.Seconds());
        scout_count = 1;
        } else {
        t.Start();
        edges_to_check -= scout_count;
        scout_count = TDStep(parent, queue);
        queue.slide_window();
        t.Stop();
    //      PrintStep("td", t.Seconds(), queue.size());
        }
    }
    #pragma omp parallel for
    for (int32_t n = 0; n < vertex_count; n++)
        if (parent[n] < -1)
        parent[n] = -1;
    return parent;
}

// // bc
// static void ParallelPrefixSum(const pvector<int32_t> &degrees, pvector<int32_t> &prefix) {
//   const size_t block_size = 1 << 20;
//   const size_t num_blocks = (degrees.size() + block_size - 1) / block_size;
//   pvector<int32_t> local_sums(num_blocks);
// #pragma omp parallel for
//   for (size_t block = 0; block < num_blocks; block++) {
//     int32_t lsum = 0;
//     size_t block_end = std::min((block + 1) * block_size, degrees.size());
//     for (size_t i = block * block_size; i < block_end; i++) lsum += degrees[i];
//     local_sums[block] = lsum;
//   }
//   pvector<int32_t> bulk_prefix(num_blocks + 1);
//   int32_t total = 0;
//   for (size_t block = 0; block < num_blocks; block++) {
//     bulk_prefix[block] = total;
//     total += local_sums[block];
//   }
//   bulk_prefix[num_blocks] = total;
// #pragma omp parallel for
//   for (size_t block = 0; block < num_blocks; block++) {
//     int32_t local_total = bulk_prefix[block];
//     size_t block_end = std::min((block + 1) * block_size, degrees.size());
//     for (size_t i = block * block_size; i < block_end; i++) {
//       prefix[i] = local_total;
//       local_total += degrees[i];
//     }
//   }
//   prefix[degrees.size()] = bulk_prefix[num_blocks];
// }

// static void ParallelLoadDegrees(pvector<int32_t> &degrees) {
// #pragma omp parallel for schedule(dynamic, 16384)
//   for (int32_t u = 0; u < global_max_vid; u++) {
//     degrees[u] = global_vectex_vec[u]->degree;
//   }
// }

// inline int32_t GetEdgeId(const pvector<int32_t> &prefix, int32_t u, int32_t local_edge_id) {
//   return prefix[u] + local_edge_id;
// }

// void PBFS(int32_t source, pvector<int32_t> &path_counts,
//     Bitmap &succ, vector<SlidingQueue<int32_t>::iterator> &depth_index,
//     SlidingQueue<int32_t> &queue, const pvector<int32_t> &prefix) {
//   pvector<int32_t> depths(global_max_vid, -1);
//   depths[source] = 0;
//   path_counts[source] = 1;
//   queue.push_back(source);
//   depth_index.push_back(queue.begin());
//   queue.slide_window();

//   #pragma omp parallel
//   {
//     int32_t depth = 0;
//     QueueBuffer<int32_t> lqueue(queue);
//     while (!queue.empty()) {
//       #pragma omp single
//       depth_index.push_back(queue.begin());
//       depth++;
//       #pragma omp for schedule(dynamic, 64)
//       for (auto q_iter = queue.begin(); q_iter < queue.end(); q_iter++) {
//         int32_t uidx = *q_iter;
//         int32_t local_edge_id = 0;

//         Vertex *u = global_vectex_vec[uidx];
//         int32_t blk_num = u->blk_list.size();
//         int32_t arr_cap = u->pma->arr_cap;
//         for (int blk_idx = (blk_num-1); blk_idx >= 0; blk_idx--) {
//           int32_t end_idx = PM_BLK_SIZE;
//           if (blk_idx == blk_num - 1 && (arr_cap % PM_BLK_SIZE != 0)) {
//             end_idx = arr_cap % PM_BLK_SIZE;
//           } 
//           int32_t *blk_ptr = u->blk_list[blk_idx].blk_ptr;
//           for (int ele_idx = (end_idx-1); ele_idx >= 0; ele_idx--) {
//             int32_t nei_v = blk_ptr[ele_idx];
//             if (nei_v >= 0) {
//               if ((depths[nei_v] == -1) && (compare_and_swap(depths[nei_v], static_cast<int32_t>(-1), depth))) {
//                 lqueue.push_back(nei_v);
//               }
//               if (depths[nei_v] == depth) {
//                 succ.set_bit_atomic(GetEdgeId(prefix, uidx, local_edge_id));
//                 fetch_and_add(path_counts[nei_v], path_counts[uidx]);
//               }
//               local_edge_id += 1;
//             }
//           }
//         }
//       }
//       lqueue.flush();
//       #pragma omp barrier
//       #pragma omp single
//       queue.slide_window();
//     }
//   }
//   depth_index.push_back(queue.begin());
// }


// pvector<ScoreT> Brandes(int32_t source, int32_t num_iters) {
//   pvector<ScoreT> scores(global_max_vid, 0);
//   pvector<int32_t> path_counts(global_max_vid);
//   Bitmap succ(global_edge_num);
//   vector<SlidingQueue<int32_t>::iterator> depth_index;
//   SlidingQueue<int32_t> queue(global_max_vid);

//   pvector<int32_t> degrees(global_max_vid);
//   ParallelLoadDegrees(degrees);
//   pvector<int32_t> prefix(degrees.size() + 1);
//   ParallelPrefixSum(degrees, prefix);

//   for (int32_t iter=0; iter < num_iters; iter++) {
//     int32_t source = source;
//     path_counts.fill(0);
//     depth_index.resize(0);
//     queue.reset();
//     succ.reset();
//     PBFS(source, path_counts, succ, depth_index, queue, prefix);
//     pvector<ScoreT> deltas(global_max_vid, 0);
//     for (int d=depth_index.size()-2; d >= 0; d--) {
//       #pragma omp parallel for schedule(dynamic, 64)
//       for (auto it = depth_index[d]; it < depth_index[d+1]; it++) {
//         int32_t uidx = *it;
//         ScoreT delta_u = 0;
//         int32_t local_edge_id = 0;

//         Vertex *u = global_vectex_vec[uidx];
//         int32_t blk_num = u->blk_list.size();
//         int32_t arr_cap = u->pma->arr_cap;
//         for (int blk_idx = (blk_num-1); blk_idx >= 0; blk_idx--) {
//           int32_t end_idx = PM_BLK_SIZE;
//           if (blk_idx == blk_num - 1 && (arr_cap % PM_BLK_SIZE != 0)) {
//             end_idx = arr_cap % PM_BLK_SIZE;
//           } 
//           int32_t *blk_ptr = u->blk_list[blk_idx].blk_ptr;
//           for (int ele_idx = (end_idx-1); ele_idx >= 0; ele_idx--) {
//             int32_t nei_v = blk_ptr[ele_idx];
//             if (nei_v >= 0) {
//               if (succ.get_bit(GetEdgeId(prefix, uidx, local_edge_id))) {
//                 delta_u += static_cast<ScoreT>(path_counts[uidx]) / static_cast<ScoreT>(path_counts[nei_v]) * (1 + deltas[nei_v]);
//               }
//               local_edge_id += 1;
//             }
//           }
//         }

//         deltas[uidx] = delta_u;
//         scores[uidx] += delta_u;
//       }
//     }
//   }
//   // normalize scores
//   ScoreT biggest_score = 0;
//   #pragma omp parallel for reduction(max : biggest_score)
//   for (int32_t n=0; n < global_max_vid; n++)
//     biggest_score = max(biggest_score, scores[n]);
//   #pragma omp parallel for
//   for (int32_t n=0; n < global_max_vid; n++)
//     scores[n] = scores[n] / biggest_score;
//   return scores;
// }


// cc
pvector<int32_t> ShiloachVishkin(void) {
    unique_lock<std::shared_timed_mutex> graph_lock(global_query_update_mutex);
    const int32_t vertex_count = QueryVertexCount();
    pvector<int32_t> comp(vertex_count);
    #pragma omp parallel for
    for (int32_t n=0; n < vertex_count; n++)
        comp[n] = n;
    bool change = true;
    int num_iter = 0;
    while (change) {
        std::atomic<bool> any_change(false);
        num_iter++;
        // note: this gives better scaleup performance
        // #pragma omp parallel for schedule(dynamic, 64)
        #pragma omp parallel for
        for (int32_t uidx=0; uidx < vertex_count; uidx++) {
            Vertex *v = global_vectex_vec[uidx];
            for_each_neighbor_readonly(v, [&](int32_t neighbor) {
                    int32_t comp_u = comp[uidx];
                    int32_t comp_v = comp[neighbor];
                    if (comp_u == comp_v) return;
                    // Hooking condition so lower component ID wins independent of direction
                    int32_t high_comp = comp_u > comp_v ? comp_u : comp_v;
                    int32_t low_comp = comp_u < comp_v ? comp_u : comp_v;
                    int32_t expected = high_comp;
                    if (compare_and_swap(comp[high_comp], expected, low_comp)) {
                        any_change.store(true, std::memory_order_relaxed);
                    }
            });
        }
        #pragma omp parallel for
        for (int32_t n=0; n < vertex_count; n++) {
            while (comp[n] != comp[comp[n]]) {
                comp[n] = comp[comp[n]];
            }
        }
        change = any_change.load(std::memory_order_relaxed);
    }
    return comp;
}

void run_query_benchmarks(int32_t source, int pagerank_iters, const string &label = "")
{
    double start_t, end_t;
    string prefix = label.empty() ? "QUERY" : "QUERY_" + label;
    const int32_t vertex_count = QueryVertexCount();
    if (vertex_count <= 0) {
        cerr << "QUERY_ERROR empty graph" << endl;
        return;
    }
    if (source < 0 || source >= vertex_count) {
        source = 1 < vertex_count ? 1 : 0;
    }

    vector<QueryPrefetchedBlock> prefetched =
        prefetch_high_degree_query_blocks(query_prefetch_block_limit);

    start_t = get_current_time();
    pvector<int32_t> bfs_parent = DOBFS(source);
    end_t = get_current_time();
    cout << prefix << "_BFS_SEC=" << (end_t - start_t) << endl;

    start_t = get_current_time();
    ScoreT *pagerank_scores = PageRankPullGS(pagerank_iters);
    end_t = get_current_time();
    cout << prefix << "_PAGERANK_SEC=" << (end_t - start_t) << endl;
    free(pagerank_scores);

    start_t = get_current_time();
    pvector<int32_t> components = ShiloachVishkin();
    end_t = get_current_time();
    cout << prefix << "_WCC_SEC=" << (end_t - start_t) << endl;

    release_query_prefetched_blocks(prefetched);
}

void test_benchmark(void)
{
    double start_t, end_t;
    global_max_vid += 1;
    for (int i = 0; i < 10; i++) {
        // cout << "PageRank start." << endl;
        // start_t = get_current_time();
        // PageRankPullGS(20);
        // end_t = get_current_time();
        // cout << "PageRank time: " << end_t - start_t << endl;

        cout << "BFS start." << endl;
        start_t = get_current_time();
        DOBFS(1);
        end_t = get_current_time();
        cout << "BFS time: " << end_t - start_t << endl;
        
    }

}

#endif 
