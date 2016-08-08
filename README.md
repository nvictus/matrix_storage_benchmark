## Matrix storage formats benchmark

Test how well each contact matrix storage format performs.

**Data Example:**

```
[peter@dbmipkedjievmbp matrix_storage_benchmark] [master]$ gzcat test/data/chrX_5KB_bins.tsv.gz | head -n 3
12      12      294.0
12      13      3.0
13      13      6.0
```

**Usage Example:**


Binned chrX example.

```
python scripts/matrix_storage.py test/data/chrX_5KB_bins.3col.sorted.tsv.gz
```

