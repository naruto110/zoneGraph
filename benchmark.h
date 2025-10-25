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
ScoreT * PageRankPullGS(int max_iters, double epsilon = 1e-4) 
{
    omp_set_num_threads(THREAD_NUM);
    const ScoreT init_score = 1.0f / global_max_vid;
    const ScoreT base_score = (1.0f - kDamp) / global_max_vid;
    ScoreT *scores;
    ScoreT *outgoing_contrib;
    scores = (ScoreT *) malloc(sizeof(ScoreT) * global_max_vid);
    outgoing_contrib = (ScoreT *) malloc(sizeof(ScoreT) * global_max_vid);

    #pragma omp parallel for
    // #pragma omp parallel for num_threads(THREAD_NUM)
    for (int32_t vidx = 0; vidx < global_max_vid; vidx++)
        outgoing_contrib[vidx] = init_score;


    for (int iter=0; iter < max_iters; iter++) {
        double error = 0;
        #pragma omp parallel for reduction(+ : error) schedule(dynamic, THREAD_NUM)
        for (int32_t vidx=0; vidx < global_max_vid; vidx++) {
        // printf("I am thread %d / %d \n", omp_get_thread_num(),omp_get_num_threads());
        // cout << "thread: " << omp_get_thread_num() << " total threads: " << omp_get_num_threads();
            ScoreT incoming_total = 0;
            Vertex *v = global_vectex_vec[vidx];
            int32_t blk_cnt = v->blk_list.size();
            int32_t *tmp_blk_ptr = nullptr;

            if (v->degree < CACHE_BLK_SIZE) {
                for (int32_t i = 0; i < v->degree; i++) {
                    incoming_total += outgoing_contrib[v->cache_buf[i]];
                }
            } else {
                if (v->blk_list.size() == 1 && v->blk_list[0].buf == nullptr) {
                    cacheBlock(v, 0);  
                }
                for (int blkidx = 0; blkidx < v->blk_list.size(); blkidx++) {
                    if (v->blk_list[blkidx].buf != nullptr) {
                        tmp_blk_ptr = v->blk_list[blkidx].buf;
                        for (int j = 0; j < PM_BLK_SIZE; j++) {
                            if (tmp_blk_ptr[j] >= 0) {
                                incoming_total += outgoing_contrib[tmp_blk_ptr[j]];
                            }
                        }
                    } else {
                        readBlockFromZone(zns.dev_fd, tmp_blk_ptr, PM_BLK_SIZE, v->blk_list[blkidx].global_wp);
                        for (int j = 0; j < PM_BLK_SIZE; j++) {
                            if (tmp_blk_ptr[j] >= 0) {
                                incoming_total += outgoing_contrib[tmp_blk_ptr[j]];
                            }
                        }
                        delete[] tmp_blk_ptr;
                        tmp_blk_ptr = nullptr;
                    }
                }
            }

            ScoreT old_score = scores[vidx];
            scores[vidx] = base_score + kDamp * incoming_total;
            error += fabs(scores[vidx] - old_score);
            outgoing_contrib[vidx] = scores[vidx] / v->degree;
            // delete[] tmp_blk_ptr;
        }
        // cout << "iter: " << iter << " total_num_access: " << total_num_access << endl;
    //  printf(" %2d    %lf\n", iter, error);
    //  if (error < epsilon)
    //    break;
    }

    return scores;
}



// BFS
int32_t BUStep(pvector<int32_t> &parent, Bitmap &front, Bitmap &next) {
    int32_t awake_count = 0;
    next.reset();
    #pragma omp parallel for reduction(+ : awake_count) schedule(dynamic, 1024)
    for (int32_t uidx=0; uidx < global_max_vid; uidx++) {
        if (parent[uidx] < 0) {
            Vertex *v = global_vectex_vec[uidx];
            int32_t blk_num = v->blk_list.size();
            int32_t arr_cap = v->pma->arr_cap;

            int32_t *tmp_blk_ptr = nullptr;

            if (v->degree < CACHE_BLK_SIZE) {
                for (int32_t i = 0; i < v->degree; i++) {
                    if (front.get_bit(v->cache_buf[i])) {
                        parent[uidx] = v->cache_buf[i];
                        awake_count++;
                        next.set_bit(uidx);
                        break;
                    }
                }
            } else {
                if (v->blk_list.size() == 1 && v->blk_list[0].buf == nullptr) {
                    cacheBlock(v, 0);  
                }
                for (int blkidx = 0; blkidx < v->blk_list.size(); blkidx++) {
                    if (v->blk_list[blkidx].buf != nullptr) {
                        tmp_blk_ptr = v->blk_list[blkidx].buf;
                        for (int j = 0; j < PM_BLK_SIZE; j++) {
                            if (tmp_blk_ptr[j] >= 0) {
                                if (front.get_bit(tmp_blk_ptr[j])) {
                                    parent[uidx] = tmp_blk_ptr[j];
                                    awake_count++;
                                    next.set_bit(uidx);
                                    break;
                                }
                            }
                        }
                    } else {
                        readBlockFromZone(zns.dev_fd, tmp_blk_ptr, PM_BLK_SIZE, v->blk_list[blkidx].global_wp);
                        for (int j = 0; j < PM_BLK_SIZE; j++) {
                            if (tmp_blk_ptr[j] >= 0) {
                                if (front.get_bit(tmp_blk_ptr[j])) {
                                    parent[uidx] = tmp_blk_ptr[j];
                                    awake_count++;
                                    next.set_bit(uidx);
                                    break;
                                }
                            }
                        }
                        delete[] tmp_blk_ptr;
                        tmp_blk_ptr = nullptr;
                    }
                }
            }
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
            int32_t blk_num = v->blk_list.size();
            int32_t arr_cap = v->pma->arr_cap;
            int32_t *tmp_blk_ptr = nullptr;

            if (v->degree < CACHE_BLK_SIZE) {
                for (int32_t i = 0; i < v->degree; i++) {
                    int32_t curr_val = parent[v->cache_buf[i]];
                    if (curr_val < 0) {
                        if (compare_and_swap(parent[v->cache_buf[i]], curr_val, uidx)) {
                            lqueue.push_back(v->cache_buf[i]);
                            scout_count += -curr_val;
                        }
                    }
                }
            } else {
                if (v->blk_list.size() == 1 && v->blk_list[0].buf == nullptr) {
                    cacheBlock(v, 0);  
                }
                for (int blkidx = 0; blkidx < v->blk_list.size(); blkidx++) {
                    if (v->blk_list[blkidx].buf != nullptr) {
                        tmp_blk_ptr = v->blk_list[blkidx].buf;
                        for (int j = 0; j < PM_BLK_SIZE; j++) {
                            if (tmp_blk_ptr[j] >= 0) {
                                if (tmp_blk_ptr[j] >= 0) {
                                    int32_t curr_val = parent[tmp_blk_ptr[j]];
                                    if (curr_val < 0) {
                                        if (compare_and_swap(parent[tmp_blk_ptr[j]], curr_val, uidx)) {
                                            lqueue.push_back(tmp_blk_ptr[j]);
                                            scout_count += -curr_val;
                                        }
                                    }
                                }
                            }
                        }
                    } else {
                        readBlockFromZone(zns.dev_fd, tmp_blk_ptr, PM_BLK_SIZE, v->blk_list[blkidx].global_wp);
                        for (int j = 0; j < PM_BLK_SIZE; j++) {
                            if (tmp_blk_ptr[j] >= 0) {
                                if (tmp_blk_ptr[j] >= 0) {
                                    int32_t curr_val = parent[tmp_blk_ptr[j]];
                                    if (curr_val < 0) {
                                        if (compare_and_swap(parent[tmp_blk_ptr[j]], curr_val, uidx)) {
                                            lqueue.push_back(tmp_blk_ptr[j]);
                                            scout_count += -curr_val;
                                        }
                                    }
                                }
                            }
                        }
                        delete[] tmp_blk_ptr;
                        tmp_blk_ptr = nullptr;
                    }
                }
            }
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
    #pragma omp parallel
    {
        QueueBuffer<int32_t> lqueue(queue);
        #pragma omp for
        for (int32_t n=0; n < global_max_vid; n++)
        if (bm.get_bit(n))
            lqueue.push_back(n);
        lqueue.flush();
    }
    queue.slide_window();
}

pvector<int32_t> InitParent(void) {
    pvector<int32_t> parent(global_max_vid);
    #pragma omp parallel for
    for (int32_t n=0; n < global_max_vid; n++)
        parent[n] = global_vectex_vec[n]->degree != 0 ? -global_vectex_vec[n]->degree : -1;
    return parent;
}

pvector<int32_t> DOBFS(int32_t source, int alpha = 15, int beta = 18) {
//  PrintStep("Source", static_cast<int32_t>(source));
    Timer t;
    t.Start();
    pvector<int32_t> parent = InitParent();
    t.Stop();
    //  PrintStep("i", t.Seconds());
    parent[source] = source;
    SlidingQueue<int32_t> queue(global_max_vid);
    queue.push_back(source);
    queue.slide_window();
    Bitmap curr(global_max_vid);
    curr.reset();
    Bitmap front(global_max_vid);
    front.reset();
    int64_t edges_to_check = global_edge_num;
    int64_t scout_count = global_vectex_vec[source]->degree;
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
                (awake_count > global_max_vid / beta));
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
    for (int32_t n = 0; n < global_max_vid; n++)
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
    pvector<int32_t> comp(global_max_vid);
    #pragma omp parallel for
    for (int32_t n=0; n < global_max_vid; n++)
        comp[n] = n;
    bool change = true;
    int num_iter = 0;
    while (change) {
        change = false;
        num_iter++;
        // note: this gives better scaleup performance
        // #pragma omp parallel for schedule(dynamic, 64)
        #pragma omp parallel for
        for (int32_t uidx=0; uidx < global_max_vid; uidx++) {
            Vertex *v = global_vectex_vec[uidx];
            int32_t blk_num = v->blk_list.size();
            int32_t arr_cap = v->pma->arr_cap;
            int32_t *tmp_blk_ptr = nullptr;

            if (v->degree < CACHE_BLK_SIZE) {
                for (int32_t i = 0; i < v->degree; i++) {
                    int32_t comp_u = comp[uidx];
                    int32_t comp_v = comp[v->cache_buf[i]];
                    if (comp_u == comp_v) continue;
                    // Hooking condition so lower component ID wins independent of direction
                    int32_t high_comp = comp_u > comp_v ? comp_u : comp_v;
                    int32_t low_comp = comp_u + (comp_v - high_comp);
                    if (high_comp == comp[high_comp]) {
                        change = true;
                        comp[high_comp] = low_comp;
                    }
                }
            } else {
                if (v->blk_list.size() == 1 && v->blk_list[0].buf == nullptr) {
                    cacheBlock(v, 0);  
                }
                for (int blkidx = 0; blkidx < v->blk_list.size(); blkidx++) {
                    if (v->blk_list[blkidx].buf != nullptr) {
                        tmp_blk_ptr = v->blk_list[blkidx].buf;
                        for (int j = 0; j < PM_BLK_SIZE; j++) {
                            if (tmp_blk_ptr[j] >= 0) {
                                int32_t comp_u = comp[uidx];
                                int32_t comp_v = comp[tmp_blk_ptr[j]];
                                if (comp_u == comp_v) continue;
                                // Hooking condition so lower component ID wins independent of direction
                                int32_t high_comp = comp_u > comp_v ? comp_u : comp_v;
                                int32_t low_comp = comp_u + (comp_v - high_comp);
                                if (high_comp == comp[high_comp]) {
                                    change = true;
                                    comp[high_comp] = low_comp;
                                }
                            }
                        }
                    } else {
                        readBlockFromZone(zns.dev_fd, tmp_blk_ptr, PM_BLK_SIZE, v->blk_list[blkidx].global_wp);
                        for (int j = 0; j < PM_BLK_SIZE; j++) {
                            if (tmp_blk_ptr[j] >= 0) {
                                int32_t comp_u = comp[uidx];
                                int32_t comp_v = comp[tmp_blk_ptr[j]];
                                if (comp_u == comp_v) continue;
                                // Hooking condition so lower component ID wins independent of direction
                                int32_t high_comp = comp_u > comp_v ? comp_u : comp_v;
                                int32_t low_comp = comp_u + (comp_v - high_comp);
                                if (high_comp == comp[high_comp]) {
                                    change = true;
                                    comp[high_comp] = low_comp;
                                }
                            }
                        }
                        delete[] tmp_blk_ptr;
                        tmp_blk_ptr = nullptr;
                    }
                }
            }
        }
        #pragma omp parallel for
        for (int32_t n=0; n < global_max_vid; n++) {
            while (comp[n] != comp[comp[n]]) {
                comp[n] = comp[comp[n]];
            }
        }
    }
    cout << "Shiloach-Vishkin took " << num_iter << " iterations" << endl;
    return comp;
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


