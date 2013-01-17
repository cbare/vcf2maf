cp foo.maf /Users/chris/Documents/work/projects/vcf2maf/soundcheck/ucsc.edu_OV.IlluminaGA-DNASeq.Level_2.128.0.0/foo.maf
cd /Users/chris/Documents/work/projects/vcf2maf/soundcheck/ucsc.edu_OV.IlluminaGA-DNASeq.Level_2.128.0.0
md5 -r * > MANIFEST.txt
cd ..
tar -czvf ucsc.edu_OV.IlluminaGA-DNASeq.Level_2.128.0.0.tar.gz ucsc.edu_OV.IlluminaGA-DNASeq.Level_2.128.0.0
./validate.sh ucsc.edu_OV.IlluminaGA-DNASeq.Level_2.128.0.0.tar.gz -bypass -noremote -usebarcode -v -centertype GSC

