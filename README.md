# ChIP-seq

All ChIP-seq data were processed according to the ENCODE3 ChIP-seq pipeline (https://www.encodeproject.org/chip-seq/), and mapped to hg19; all data passed ENCODE quality standards. ChIP-seq peaks were called using MACS2 (https://github.com/taoliu/MACS), followed by identifying common peaks between duplicates using IDR (https://github.com/nboley/idr).

![image](https://github.com/swniw/ChIP-seq/assets/120678930/eb54b0ea-339e-4443-89df-de5701de8312)

Reference: 
- [ENCODE Uniform Processing Pipelines](https://www.encodeproject.org/chip-seq/transcription_factor/)
- [W. Ni, A. A. Perez, S. Schreiner, C. M. Nicolet, and P. J. Farnham. Characterization of the ZFX family of transcription factors that bind downstream of the start site of CpG island promoters. Nucleic Acids Res, 2020. doi: 10.1093/nar/gkaa384.](https://academic.oup.com/nar/article/48/11/5986/5837054?login=false)
