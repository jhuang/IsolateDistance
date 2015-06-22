#!/usr/bin/env ruby

require "pp"
require "open3"
require "set"
# TODO: add freebayes (too much SNPs) or GATK
# require "options"


if ARGV.length == 0 then
  STDERR.puts "Usage: #{$0} *.fastq.gz"
  STDERR.puts "Example: #{$0} *_R1_*.fastq.gz"
  STDERR.puts "Example: #{$0} ./dataset_test/1-VA1290-h-rohde_S7_L001_R1_001_l500000.fastq.gz ./dataset_test/2-VA4851-h-rohde_S8_L001_R1_001_l500000.fastq.gz"
  STDERR.puts "Example: #{$0} 1-VA1290-h-rohde_S7_L001_R1_001.fastq.gz 2-VA4851-h-rohde_S8_L001_R1_001.fastq.gz \
  3-VA56187-h-rohde_S9_L001_R1_001.fastq.gz 4-VA1620-h-rohde_S10_L001_R1_001.fastq.gz \
  5-VA1664-h-rohde_S11_L001_R1_001.fastq.gz 6-VA2382-h-rohde_S12_L001_R1_001.fastq.gz \
  7-VA4392-h-rohde_S13_L001_R1_001.fastq.gz 8-VA8509-h-rohde_S14_L001_R1_001.fastq.gz \
  9-VA9248-h-rohde_S15_L001_R1_001.fastq.gz"
  STDERR.puts "Example: #{$0} ./dataset1/ERR101899_1.fastq.gz ./dataset1/ERR101900_1.fastq.gz \
  ./dataset1/ERR103394_1.fastq.gz ./dataset1/ERR103395_1.fastq.gz ./dataset1/ERR103396_1.fastq.gz ./dataset1/ERR103397_1.fastq.gz \
  ./dataset1/ERR103398_1.fastq.gz ./dataset1/ERR103400_1.fastq.gz ./dataset1/ERR103401_1.fastq.gz ./dataset1/ERR103402_1.fastq.gz \
  ./dataset1/ERR103403_1.fastq.gz ./dataset1/ERR103404_1.fastq.gz ./dataset1/ERR103405_1.fastq.gz ./dataset1/ERR159680_1.fastq.gz"
  STDERR.puts "Example: #{$0} ./dataset2/ERR435878_1.fastq.gz ./dataset2/ERR436034_1.fastq.gz \
  ./dataset2/ERR436035_1.fastq.gz ./dataset2/ERR436036_1.fastq.gz ./dataset2/ERR440528_1.fastq.gz ./dataset2/ERR440529_1.fastq.gz \
  ./dataset2/ERR454986_1.fastq.gz ./dataset2/ERR458146_1.fastq.gz ./dataset2/ERR458147_1.fastq.gz ./dataset2/ERR468930_1.fastq.gz \
  ./dataset2/ERR469300_1.fastq.gz ./dataset2/ERR469301_1.fastq.gz ./dataset2/ERR469619_1.fastq.gz ./dataset2/ERR469620_1.fastq.gz \
  ./dataset2/ERR469621_1.fastq.gz ./dataset2/ERR469622_1.fastq.gz ./dataset2/ERR469623_1.fastq.gz ./dataset2/ERR469624_1.fastq.gz \
  ./dataset2/ERR469625_1.fastq.gz ./dataset2/ERR469626_1.fastq.gz ./dataset2/ERR469627_1.fastq.gz"
  exit(-1)
end

#Path to Bowtie2
SPADES_PATH="/home/jhuang/Tools/SPAdes-3.5.0-Linux"
SPADES="/home/jhuang/Tools/SPAdes-3.5.0-Linux/bin/spades.py"
BOWTIE2_BUILD="/home/jhuang/Tools/bowtie2-2.2.5/bowtie2-build"
BOWTIE2="/home/jhuang/Tools/bowtie2-2.2.5/bowtie2"
BWA="/usr/local/bin/bwa";
VCFFILTER="/home/jhuang/Tools/freebayes_git/freebayes/vcflib/bin/vcffilter"
FREEBAYES="/usr/bin/freebayes"
CRESEQDIC="picard-tools CreateSequenceDictionary"
GATK="/home/jhuang/Tools/GenomeAnalysisTK-3.4-0/GenomeAnalysisTK.jar"
GATK2="/home/jhuang/Tools/GenomeAnalysisTK-2.8-1/GenomeAnalysisTK.jar"
SAMTOOLS="/usr/bin/samtools"
#ReduceReads




#### Step 0: assemblies and mappings ####
assemblies = Array.new
read_1s = ARGV  # with relative path
read_2s = []  # with relative path
rel_paths = []
sample_names = []  # without relative path
sample_name_id = {}  # without relative path
read_1s.each_with_index do |read_1, idx|
  if read_1.include?("_R1_") then
    read_2 = read_1.sub("_R1_","_R2_")
  elsif read_1.include?("_1") then
    read_2 = read_1.sub("_1","_2")
  else
    puts "The pattern of fastq.gz files is not recognized."
    exit(-1)
  end
  read_2s << read_2
  
  assembly = "#{read_1.chomp(".fastq.gz")}_assembly/contigs.fasta"
  if !File.exist?(assembly) then
    puts "running spades.py --careful --threads 15 -1 #{read_1} -2 #{read_2} -o #{read_1.chomp(".fastq.gz")}_assembly"
    `#{SPADES} --careful --threads 15 -1 #{read_1} -2 #{read_2} -o #{read_1.chomp(".fastq.gz")}_assembly`
    `samtools faidx #{assembly}`
  end
  assemblies << assembly
  
  if read_1.split("/").length >= 2 then
    rel_paths << read_1.split("/")[0..-2].join("/")
  else
    rel_paths << ""
  end
  sample_name = read_1.split("/")[-1].chomp(".fastq.gz")
  sample_name_id[sample_name] = idx 
  sample_names << sample_name
end
pp rel_paths
if rel_paths.uniq.length != 1 then
  puts "The relative path of input files have to be the same!"
  exit(-1)
end
rel_path = rel_paths.uniq[0]
puts "rel_path='#{rel_path}'"


#### 1. direct generate the merged_bams ####
# check if a readid exists in all mappings

# a read from sample j, it should all mapped on different assemblies[i]
# check if it exists on all mappings
0.upto(sample_names.length-1) do |i|
  # TODO: make more strict check in the following line
#   if File.exist?("#{rel_path}/#{sample_names[i]}-#{sample_names[0]}_filtered.bam") then
#     break
#   end

  # prepare read1 filepath and read2 filepath: reads from sample i
  read1_file = "#{assemblies[i].sub("contigs.fasta", "corrected")}/#{read_1s[i].sub("#{rel_path}/", "").chomp(".gz")}.00.0_0.cor.fastq.gz"
  read2_file = "#{assemblies[i].sub("contigs.fasta", "corrected")}/#{read_2s[i].sub("#{rel_path}/", "").chomp(".gz")}.00.0_0.cor.fastq.gz"
  
  h = {}
  assemblies.each_with_index do |assembly, j| 
    ## Bowtie2 ##
    # TODO: using fork to process data simultaneously
    if j==0 then
      puts "#{BOWTIE2_BUILD} #{assemblies[j]} #{assemblies[j]}_index"
      `#{BOWTIE2_BUILD} #{assemblies[j]} #{assemblies[j]}_index`
      puts "#{BOWTIE2} -x #{assemblies[j]}_index -1 #{read1_file} -2 #{read2_file} --threads 15 --rg ID:#{sample_names[i]} --rg SM:#{sample_names[i]} --rg PL:illumina --rg LB:standard --very-sensitive"
      #`#{BOWTIE2} -x #{assemblies[j]}_index -1 #{read1_file} -2 #{read2_file} --threads 15 --rg ID:#{sample_names[i]} --rg SM:#{sample_names[i]} --rg PL:illumina --rg LB:standard --very-sensitive | samtools view -bS -f0x2 - | samtools sort -n - #{rel_path}/#{sample_names[i]}-#{sample_names[j]}_sorted`
      `#{BOWTIE2} -x #{assemblies[j]}_index -1 #{read1_file} -2 #{read2_file} --threads 15 --rg ID:#{sample_names[i]} --rg SM:#{sample_names[i]} --rg PL:illumina --rg LB:standard --very-sensitive | samtools view -bS - | samtools sort -n - #{rel_path}/#{sample_names[i]}-#{sample_names[j]}_sbyn`
      `samtools view -hf 2 #{rel_path}/#{sample_names[i]}-#{sample_names[j]}_sbyn.bam | grep -v "XS:i:" | ./foo.py | samtools view - -Sbo #{rel_path}/#{sample_names[i]}-#{sample_names[j]}_sorted.bam`  ## Note that the file only contains Reads-Pairs Aligned Concordantly Exactly 1 Time --> TODO: correct this point.
      f = IO.popen("bamToBed -i #{rel_path}/#{sample_names[i]}-#{sample_names[j]}_sorted.bam -bedpe | cut -f7")
      f.each_line do |line|
	line = line.chomp
	h[line] = 1
      end
    else
      puts "#{BOWTIE2_BUILD} #{assemblies[j]} #{assemblies[j]}_index"
      `#{BOWTIE2_BUILD} #{assemblies[j]} #{assemblies[j]}_index`
      puts "#{BOWTIE2} -x #{assemblies[j]}_index -1 #{read1_file} -2 #{read2_file} --threads 15 --rg ID:#{sample_names[i]} --rg SM:#{sample_names[i]} --rg PL:illumina --rg LB:standard --very-sensitive"
      #`#{BOWTIE2} -x #{assemblies[j]}_index -1 #{read1_file} -2 #{read2_file} --threads 15 --rg ID:#{sample_names[i]} --rg SM:#{sample_names[i]} --rg PL:illumina --rg LB:standard --very-sensitive | samtools view -bS -f0x2 - | samtools sort -n - #{rel_path}/#{sample_names[i]}-#{sample_names[j]}_sorted`
      `#{BOWTIE2} -x #{assemblies[j]}_index -1 #{read1_file} -2 #{read2_file} --threads 15 --rg ID:#{sample_names[i]} --rg SM:#{sample_names[i]} --rg PL:illumina --rg LB:standard --very-sensitive | samtools view -bS - | samtools sort -n - #{rel_path}/#{sample_names[i]}-#{sample_names[j]}_sbyn`
      `samtools view -hf 2 #{rel_path}/#{sample_names[i]}-#{sample_names[j]}_sbyn.bam | grep -v "XS:i:" | ./foo.py | samtools view - -Sbo #{rel_path}/#{sample_names[i]}-#{sample_names[j]}_sorted.bam`  ## Note that the file only contains Reads-Pairs Aligned Concordantly Exactly 1 Time.
      f = IO.popen("bamToBed -i #{rel_path}/#{sample_names[i]}-#{sample_names[j]}_sorted.bam -bedpe | cut -f7")
      f.each_line do |line|
	line = line.chomp
	if h[line] && h[line]==j then
	  h[line] += 1
	elsif h[line] then
	  h.delete(line)
	end
      end    
    end
    puts "h.size=#{h.size}"
  end
  #pp h

  
  assemblies.each_with_index do |assembly, j|
    ## from {sample_names[i]}-#{sample_names[j]}_sorted.bam to #{sample_names[i]}-#{sample_names[j]}_filtered.bam
    f_filtered_sam = File.open("#{rel_path}/#{sample_names[i]}-#{sample_names[j]}_filtered.sam", "w")
    puts "IO.popen('samtools view -h #{rel_path}/#{sample_names[i]}-#{sample_names[j]}_sorted.bam')"
    f = IO.popen("samtools view -h #{rel_path}/#{sample_names[i]}-#{sample_names[j]}_sorted.bam")
    #idx = 0
    f.each_line do |line|
      #if idx%100000==0 then
      #  puts idx
      #end
      if line.start_with?("@") then
	f_filtered_sam.puts(line)
      else
# 	# inefficient solution
# 	token = line.split(/\t/)
# 	readid = token.at(0)
# 	if filtered_reads.include?(readid.chomp) then
# 	    f_filtered_sam.print(line)
# 	end
	if (m = /(\S+)\s+/.match(line)) && h[m[1].chomp] then
	  f_filtered_sam.puts(line)
	else
# 	  puts h[m[1].chomp]
# 	  puts m[1].chomp
	end
      end
      #idx += 1
    end
    f.close
    f_filtered_sam.close
    #samtools view -bS 1-VA1290-h-rohde_S7_L001_R1_001-1-VA1290-h-rohde_S7_L001_R1_001_filtered.sam > 1-VA1290-h-rohde_S7_L001_R1_001-1-VA1290-h-rohde_S7_L001_R1_001_filtered.bam
    `samtools view -bS #{rel_path}/#{sample_names[i]}-#{sample_names[j]}_filtered.sam > #{rel_path}/#{sample_names[i]}-#{sample_names[j]}_filtered.bam`
    `rm #{rel_path}/#{sample_names[i]}-#{sample_names[j]}_filtered.sam`
  end
end


merged_bams = []
assemblies.each_with_index do |assembly, j| 
#   if !File.exist?("#{rel_path}/#{sample_names[j]}_merged_sorted.bam") then
    # fix 'samtools merge' ignoring @RG
    `samtools view -H #{rel_path}/#{sample_names[0]}-#{sample_names[j]}_sorted.bam | sed '/@RG/d' - > #{rel_path}/#{sample_names[j]}_header.sam`  
    filtered_bams = []
    
    0.upto(sample_names.length-1) do |i|  # read from different samples ==> marked with RG
      filtered_bams << "#{rel_path}/#{sample_names[i]}-#{sample_names[j]}_filtered.bam"
      `echo "@RG\tID:#{sample_names[i]}\tSM:#{sample_names[i]}\tPL:illumina\tLB:standard" >> #{rel_path}/#{sample_names[j]}_header.sam`
    end
    #`samtools merge -rh #{sample_names[j]}_header.sam #{sample_names[j]}_merged.bam #{filtered_bams.join(" ")}`
    if File.exist?("#{rel_path}/#{sample_names[j]}_merged.bam") then
      `rm #{rel_path}/#{sample_names[j]}_merged.bam`
    end
    `samtools merge -h #{rel_path}/#{sample_names[j]}_header.sam #{rel_path}/#{sample_names[j]}_merged.bam #{filtered_bams.join(" ")}`
    `samtools sort #{rel_path}/#{sample_names[j]}_merged.bam #{rel_path}/#{sample_names[j]}_merged_sorted`
    `samtools index #{rel_path}/#{sample_names[j]}_merged_sorted.bam`
    `rm #{rel_path}/#{sample_names[j]}_header.sam`
    `rm #{rel_path}/#{sample_names[j]}_merged.bam`
#   end  
  merged_bams << "#{rel_path}/#{sample_names[j]}_merged_sorted.bam"  # merged_bams are with relative path
end


#### Step 4: SNP calling and filtering ####
merged_bams.each_with_index do |bam_filename, j|
#   if !File.exist?("#{rel_path}/#{sample_names.at(j)}_passed.vcf") then  
    if !File.exist?(assemblies.at(j).sub(".fasta", ".dict")) then
      `#{CRESEQDIC} R=#{assemblies.at(j)} O=#{assemblies.at(j).sub(".fasta", ".dict")}`
    end
    #### using FREEBAYES ####
#     puts "#{FREEBAYES} --fasta-reference #{assemblies.at(j)} -b #{bam_filename} -p 1 -i -X -u | #{VCFFILTER} -f 'QUAL > 30' > #{rel_path}/#{sample_names.at(j)}_raw.vcf"
#     `#{FREEBAYES} --fasta-reference #{assemblies.at(j)} -b #{bam_filename} -p 1 -i -X -u | #{VCFFILTER} -f 'QUAL > 30' > #{rel_path}/#{sample_names.at(j)}_raw.vcf`

    ### using GATK ####
      puts "java -jar #{GATK} -glm SNP -R #{sample_names.at(j)}_assembly/contigs.fasta -T UnifiedGenotyper -I #{bam_filename} -o #{sample_names.at(j)}_raw.vcf -ploidy 1 -L #{sample_names.at(j)}.bed"
      `java -jar #{GATK} -glm SNP -R #{sample_names.at(j)}_assembly/contigs.fasta -T UnifiedGenotyper -I #{bam_filename} -o #{sample_names.at(j)}_raw.vcf -ploidy 1 -L #{sample_names.at(j)}.bed`
      puts "java -jar #{GATK} -glm SNP -R #{sample_names.at(j)}_assembly/contigs.fasta -T UnifiedGenotyper -I #{bam_filename} -o #{sample_names.at(j)}_raw.vcf -ploidy 1"
      `java -jar #{GATK} -glm SNP -R #{sample_names.at(j)}_assembly/contigs.fasta -T UnifiedGenotyper -I #{bam_filename} -o #{sample_names.at(j)}_raw.vcf -ploidy 1`

  
    ## EVALUATE the depth of coverage:
    #   puts "java -jar #{GATK} -T DepthOfCoverage -mmq 30 -R #{sample_names.at(j)}_assembly/contigs.fasta -omitBaseOutput -o #{sample_names.at(j)}.sample_summary -I #{sample_names.at(j)}_sorted.bam"
    #   `java -jar #{GATK} -T DepthOfCoverage -mmq 30 -R #{sample_names.at(j)}_assembly/contigs.fasta -omitBaseOutput -o #{sample_names.at(j)}.sample_summary -I #{sample_names.at(j)}_sorted.bam`
    #   avg_depth=`awk /^Total/'{printf $3}' #{sample_names.at(j)}.sample_summary`
    #   
    #     --filterExpression "MLEAF < 0.95" --filterName "AFFilter"
    #     --filterExpression "QUAL < 30.0 || DP < (#{avg_depth}/2)"  --filterName "StandardFilters"`
    #   
    #     --filterExpression "MQ0 >= 4 && '((MQ0 / (1.0 * DP)) > 0.1)'" --filterName "HARD_TO_VALIDATE"

    
    ## QUAL is standards
# NODE_2_length_252035_cov_20.0204_ID_3	244727	.	G	A	14511.50	PASS	
# AB=0;ABP=0;AC=10;AF=0.714286;AN=14;AO=458;BVAR;CIGAR=1X;
# DP=643;DPB=643;DPRA=0.99027;EPP=23.663;EPPR=3.11594;GTI=0;LEN=1;
# MEANALT=1;MQM=40.5459;MQMR=38.6162;NS=14;NUMALT=1;ODDS=147.527;
# PAIRED=1;PAIREDR=1;PAO=0;PQA=0;PQR=0;PRO=0;
# QA=16497;QR=6318;RO=185;RPP=90.7036;RPPR=38.5168;
# RUN=1;SAF=317;SAP=149.874;SAR=141;SRF=122;SRP=43.8692;SRR=63;TYPE=snp;technology.illumina=1
  
    avg_depth=10
    `java -jar #{GATK} -T VariantFiltration -R #{assemblies[j]} -o #{rel_path}/#{sample_names.at(j)}_filtered.vcf -V #{rel_path}/#{sample_names.at(j)}_raw.vcf \
    --clusterSize 3 --clusterWindowSize 10 \
    --filterExpression "QUAL < 30.0"  --filterName "LQFilter" \
    --filterExpression "QD < 10.0" --filterName "QDFilter"  \
    --filterExpression "MQ < 30.0" --filterName "MQFilter"  \
    --filterExpression "FS > 10.0" --filterName "FSFilter"  \
    --filterExpression "HaplotypeScore > 20.0" --filterName "HaplotypeScoreFilter"`
    header_line_no=`grep -n "#CHROM" #{rel_path}/#{sample_names.at(j)}_filtered.vcf | cut -d':' -f 1`
    puts "head -n #{header_line_no.chomp} #{rel_path}/#{sample_names.at(j)}_filtered.vcf > #{rel_path}/#{sample_names.at(j)}.snp_head"
    `head -n #{header_line_no.chomp} #{rel_path}/#{sample_names.at(j)}_filtered.vcf > #{rel_path}/#{sample_names.at(j)}.snp_head`
    `cat #{rel_path}/#{sample_names.at(j)}_filtered.vcf | grep PASS | cat #{rel_path}/#{sample_names.at(j)}.snp_head - > #{rel_path}/#{sample_names.at(j)}_passed.vcf`
#   end
end



#### Step 5: calculating matrix ####
vcf_pos_offset = 9
width = sample_names.length
height = sample_names.length

vcfs = merged_bams.map { |x|  # vcfs are with relative path
  x.sub("_merged_sorted.bam", "_passed.vcf")
}
#ref_idx = 0
dist_matrices = []
dist_matrix_sum = Array.new(width){Array.new(height, 0)}
vcfs.each_with_index do |vcf|
  corrected_vcf = File.new("#{vcf.sub(".vcf", "_filtered.vcf")}", "w")
  #if !File.exist?(vcf) then
  f = File.new(vcf)
  dist_matrix = Array.new(width){Array.new(height, 0)}
  f.each_line do |line|
    if !line.start_with?("#") then
      token = line.split(/\t/)
      
      # delete the line that contains "."
      # delete the line having the same variant calling
      is_valid = true
      
#       if token.at(4).include?(",") then
# 	is_valid = false  
#       else
      0.upto(width-1) do |i|
	if token.at(vcf_pos_offset + i).chomp == "." then  # .split(":")[0]
	  is_valid = false
	  break
	end

# 	ad_arr = token.at(vcf_pos_offset + i).split(":")[1].split(",")
# 	if ad_arr.length != 2 then
# 	  is_valid = false
# 	  next
# 	else
# 	  if (ad_arr[0].to_f+ad_arr[1].to_f) != 0 && ad_arr[0].to_f/(ad_arr[0].to_f+ad_arr[1].to_f) < 0.8 && ad_arr[0].to_f/(ad_arr[0].to_f+ad_arr[1].to_f) > 0.2 then
# 	    ## pp ad_arr
# 	    is_valid = false
# 	    next
# 	  end
# 	end
	token_samp = token.at(vcf_pos_offset + i).split(":")
	dp = token_samp[1].to_i
	ro = token_samp[2].to_i
	ao = token_samp[4].to_i
	rf = ro/dp
	if (dp) <= 10 && rf < 0.8 && rf > 0.2 then
	  is_valid = false
	  break
	end
      end

      if is_valid then
	corrected_vcf.puts(line)
	0.upto(width-1) do |i|
	  0.upto(width-1) do |j|
	    if token.at(vcf_pos_offset + i).split(":")[0] != token.at(vcf_pos_offset + j).split(":")[0] then
	      dist_matrix_sum[i][j] += 1 
	      dist_matrix[i][j] += 1
	    end
	  end
	end
      end
    else
      corrected_vcf.puts(line)
    end
  end
  f.close
  corrected_vcf.close
  #ref_idx += 1
  dist_matrices << dist_matrix
end
pp dist_matrices


# 是有可能的，因为　Assembly 有长有短．如果一个片段　存在于　一个assembly,而不存在另一个assembly,那么Variant calling 的SNP数目就会不一样．
dist_matrix_arv = Array.new(width){Array.new(height, 0)}
0.upto(width-1) do |i|
  0.upto(width-1) do |j|
    dist_matrix_arv[i][j] = dist_matrix_sum[i][j] / width
  end
end
pp dist_matrix_arv


# Standard deviation
matrix_quaddiff_sum = Array.new(width){Array.new(height, 0)}
dist_matrices.each do |m|
  0.upto(width-1) do |i|
    0.upto(width-1) do |j|
      if (m[i][j] != 0) then
        matrix_quaddiff_sum[i][j] += (m[i][j] - dist_matrix_arv[i][j])**2
      end
    end
  end
end
sd = Array.new(width){Array.new(height, 0)}
0.upto(width-1) do |i|
  0.upto(width-1) do |j|
    sd[i][j] = (Math.sqrt(matrix_quaddiff_sum[i][j] / width)).floor
  end
end
pp sd


#### Step 6: generating distance matrix in the Phylip format and newick format using NJ ####
short_sample_names = sample_names.map {|name| name.split("-")[1]}
f = File.new("#{short_sample_names.join("-")}.dist", "w")
f.puts(" #{width}")
0.upto(width-1) do |i|
  f.print "#{sample_names[i]} "
  0.upto(width-1) do |j|
    f.print "#{dist_matrix_arv[i][j]} "
  end
  f.print "\n"
end
# #http://www.trex.uqam.ca/index.php?action=trex#methodsChoice
# # Using SplitsTree4