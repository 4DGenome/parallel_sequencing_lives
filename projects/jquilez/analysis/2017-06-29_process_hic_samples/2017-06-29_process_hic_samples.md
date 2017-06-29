# 2017-06-29_process_hic_samples
----------------------------------------------------------------------------------------------------

**objective: process Hi-C samples**

**paths are relative to /users/GR/mb/jquilez/projects/parallel_sequencing_lives**



## [2017-06-29] Run FastQC
----------------------------------------------------------------------------------------------------

```
# variables
#samples="dc3a1e069_51720e9cf b1913e6c1_51720e9cf dc3a1e069_ec92aa0bb b7fa2d8db_bfac48760"
#process=quality_control
#project=jquilez
#analysis=2017-06-28_process_hic_samples
#submit_to_cluster=no
#integrate_metadata="yes"
scripts/utils/quality_control.sh
```



## [2017-06-29] Run Hi-C pipeline
----------------------------------------------------------------------------------------------------

Configuration file used:
```


