# Troubleshooting data downloading

Occasionally, a user may not be able to connect to cistrome.org from their institutional server due to some security measure. To circumvent this, one can manually install the data required to run LISA. 

First, on the server, use the command below to issue the download URL for the required dataset.

*server*
```bash
$ lisa download --url hg38
http://cistrome.org/~alynch/data/lisa_data/hg38_2.1.tar.gz
```

Copy this URL, then on your local machine, download LISA's required data from cistrome.org (this command fetches the human genome (hg38) data for all LISA versions 2.1.x).

*local*
```bash
$ wget http://cistrome.org/~alynch/data/lisa_data/hg38_2.1.tar.gz
```

Now, transfer the downloaded LISA data from your local machine to your working directory on the server:

*local*
```bash
$ scp ./hg38_2.1.tar.gz {user}@{server}:/{PATH_TO_DIR}/
```

Once the data is has transferred, the last step is to unpack the data in the server package and install it to LISA's package directory:

*server*
```bash
(lisa_env) $ lisa unpack hg38_2.1.tar.gz --remove
```

The LISA site package folder should now contain a directory called ```data``` with the structure:
```
data/
├── hg38
│   ├── ChIP-seq_binding.npz
│   ├── Motifs_binding.npz
│   ├── RP_map.npz
│   ├── bins.bed
│   ├── gene_locs.txt
│   ├── genes.tsv
│   ├── hg38_basic_2.1_RP_map.npz
│   ├── hg38_version.txt
│   ├── lisa_data_hg38_reads_normalized.h5
│   └── metadata
│       ├── lisa_meta.tsv
│       └── motifs_meta.tsv
```