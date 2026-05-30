# Graph Experiment Results

All times are seconds. `-` means skipped or unavailable.


## Graph Maintenance / Load


| Dataset | zoneGraph | LiveGraph | LSMGraph | MBFGraph |
| --- | ---: | ---: | ---: | ---: |
| com-DBLP | 0.133 | 1.492 | 1.600 | 6.222 |
| com-Youtube | 0.434 | 5.034 | 4.600 | 11.223 |
| twitch-gamers | 1.263 | 7.084 | 12.000 | 26.208 |
| as-Skitter | 1.771 | 19.994 | 18.000 | 41.046 |
| soc-Pokec | 4.508 | 166.939 | 40.000 | 85.538 |
| soc-LiveJournal1 | 10.860 | - | 140.000 | 283.934 |


## zoneGraph Query: Before Final Dump


| Dataset | BFS | PageRank | WCC |
| --- | ---: | ---: | ---: |
| com-DBLP | 0.023 | 0.072 | 0.026 |
| com-Youtube | 0.043 | 0.228 | 0.062 |
| twitch-gamers | 0.015 | 0.165 | 0.086 |
| as-Skitter | 0.072 | 0.449 | 0.092 |
| soc-Pokec | 0.056 | 0.937 | 0.161 |
| soc-LiveJournal1 | 0.161 | 2.143 | 0.494 |


## zoneGraph Query: After Final Dump


| Dataset | BFS | PageRank | WCC |
| --- | ---: | ---: | ---: |
| com-DBLP | 0.054 | 0.136 | 0.053 |
| com-Youtube | 0.162 | 0.560 | 0.167 |
| twitch-gamers | 0.389 | 1.677 | 0.273 |
| as-Skitter | 0.537 | 1.967 | 0.325 |
| soc-Pokec | 0.949 | 4.697 | 0.514 |
| soc-LiveJournal1 | 2.042 | 10.647 | 1.421 |


## LiveGraph Query


| Dataset | BFS | PageRank | WCC |
| --- | ---: | ---: | ---: |
| com-DBLP | 0.081 | 0.340 | 0.053 |
| com-Youtube | 0.152 | 0.899 | 0.135 |
| twitch-gamers | 0.147 | 1.207 | 0.095 |
| as-Skitter | 0.153 | 2.376 | 0.638 |
| soc-Pokec | 0.543 | 6.333 | 1.409 |
| soc-LiveJournal1 | - | - | - |


## LSMGraph Query


| Dataset | BFS no-dedup | BFS dedup | PageRank no-dedup | PageRank dedup | WCC no-dedup | WCC dedup |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| com-DBLP | 0.191 | 0.181 | 0.814 | 0.899 | 0.183 | 0.206 |
| com-Youtube | 0.306 | 0.346 | 3.161 | 4.630 | 0.441 | 0.490 |
| twitch-gamers | 0.148 | 0.190 | 2.535 | 4.039 | 0.213 | 0.313 |
| as-Skitter | 0.494 | 0.614 | 5.750 | 10.816 | 0.531 | 0.766 |
| soc-Pokec | 0.412 | 0.523 | 8.875 | 12.728 | 0.653 | 0.830 |
| soc-LiveJournal1 | 3.467 | 4.017 | 50.062 | 58.420 | 3.813 | 4.410 |


## MBFGraph Query


| Dataset | BFS no-dedup | BFS dedup | PageRank no-dedup | PageRank dedup | WCC no-dedup | WCC dedup |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| com-DBLP | 13.091 | 19.236 | 9.634 | 14.221 | 13.745 | 21.160 |
| com-Youtube | 13.003 | 29.413 | 12.989 | 30.712 | 13.039 | 31.121 |
| twitch-gamers | 6.881 | 34.356 | 4.701 | 17.677 | 7.370 | 37.821 |
| as-Skitter | 19.558 | 103.510 | 17.666 | 110.097 | 19.352 | 126.143 |
| soc-Pokec | 9.991 | 119.289 | 7.820 | 118.724 | 10.167 | 147.046 |
| soc-LiveJournal1 | 22.965 | 379.842 | 18.340 | 180.366 | 24.192 | 357.887 |
