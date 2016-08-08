## Matrix storage formats benchmark

Test how well each contact matrix storage format performs.

**Data Example:**

```
[peter@dbmipkedjievmbp matrix_storage_benchmark] [master]$ gzcat test/data/chrX_5KB_bins.tsv.gz | head -n 3
12      12      294.0
12      13      3.0
13      13      6.0
```


**Binned chrX example.**

```
$ python scripts/cooler_binned.py test/data/chrX_5KB_bins.3col.sorted.tsv.gz
```

**HIC003 example.**

```
$ wget -P test/data ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1551nnn/GSM1551552/suppl/GSM1551552_HIC003_merged_nodups.txt.gz
```

To generate a sorted and indexed contacts file, fill in the `COOLER_PATH` variable of `preprocess_pairs.sh` and run

```
$ python scripts/preprocess_pairs.sh test/data/GSM1551552_HIC003_merged_nodups.txt.gz
```

Then

```
$ python scripts/cooler_pairs.py test/data/GM12878-MboI-contacts.txt.gz
```


