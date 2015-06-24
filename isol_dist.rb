#!/usr/bin/env ruby

require "pp"
require "open3"
require "set"
# TODO: add freebayes (too much SNPs) or GATK
require "./options_isol_dist.rb"


# if ARGV.length == 0 then
#   STDERR.puts "Usage: #{$0} *.fastq.gz"
#   STDERR.puts "Example: #{$0} *_R1_*.fastq.gz"
#   STDERR.puts "Example: #{$0} ./dataset_test/1-VA1290-h-rohde_S7_L001_R1_001_l500000.fastq.gz ./dataset_test/2-VA4851-h-rohde_S8_L001_R1_001_l500000.fastq.gz"
#   STDERR.puts "Example: #{$0} 1-VA1290-h-rohde_S7_L001_R1_001.fastq.gz 2-VA4851-h-rohde_S8_L001_R1_001.fastq.gz \
#   3-VA56187-h-rohde_S9_L001_R1_001.fastq.gz 4-VA1620-h-rohde_S10_L001_R1_001.fastq.gz \
#   5-VA1664-h-rohde_S11_L001_R1_001.fastq.gz 6-VA2382-h-rohde_S12_L001_R1_001.fastq.gz \
#   7-VA4392-h-rohde_S13_L001_R1_001.fastq.gz 8-VA8509-h-rohde_S14_L001_R1_001.fastq.gz \
#   9-VA9248-h-rohde_S15_L001_R1_001.fastq.gz"
#   STDERR.puts "Example: #{$0} ./dataset1/ERR101899_1.fastq.gz ./dataset1/ERR101900_1.fastq.gz \
#   ./dataset1/ERR103394_1.fastq.gz ./dataset1/ERR103395_1.fastq.gz ./dataset1/ERR103396_1.fastq.gz ./dataset1/ERR103397_1.fastq.gz \
#   ./dataset1/ERR103398_1.fastq.gz ./dataset1/ERR103400_1.fastq.gz ./dataset1/ERR103401_1.fastq.gz ./dataset1/ERR103402_1.fastq.gz \
#   ./dataset1/ERR103403_1.fastq.gz ./dataset1/ERR103404_1.fastq.gz ./dataset1/ERR103405_1.fastq.gz ./dataset1/ERR159680_1.fastq.gz"
#   STDERR.puts "Example: #{$0} ./dataset2/ERR435878_1.fastq.gz ./dataset2/ERR436034_1.fastq.gz \
#   ./dataset2/ERR436035_1.fastq.gz ./dataset2/ERR436036_1.fastq.gz ./dataset2/ERR440528_1.fastq.gz ./dataset2/ERR440529_1.fastq.gz \
#   ./dataset2/ERR454986_1.fastq.gz ./dataset2/ERR458146_1.fastq.gz ./dataset2/ERR458147_1.fastq.gz ./dataset2/ERR468930_1.fastq.gz \
#   ./dataset2/ERR469300_1.fastq.gz ./dataset2/ERR469301_1.fastq.gz ./dataset2/ERR469619_1.fastq.gz ./dataset2/ERR469620_1.fastq.gz \
#   ./dataset2/ERR469621_1.fastq.gz ./dataset2/ERR469622_1.fastq.gz ./dataset2/ERR469623_1.fastq.gz ./dataset2/ERR469624_1.fastq.gz \
#   ./dataset2/ERR469625_1.fastq.gz ./dataset2/ERR469626_1.fastq.gz ./dataset2/ERR469627_1.fastq.gz"
#   exit(-1)
# end

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
MUGSY="mugsy"
#ReduceReads



options = IsolDistOptparser.parse(ARGV)
pp options
inputfiles = options.inputfiles
outputdir = options.outputdir
if (inputfiles.empty? || outputdir.empty?)
  puts ""
  puts "--- Please specify input files and output dir ---"
  puts ""
  puts "Please type #{$0} -h for help in usage."
  puts ""
  #puts "Could not detect the input file, please use option -i to read input. See -h for more possible usages."
  exit 1
end

inputfiles = options.inputfiles
outputdir = options.outputdir
steps = options.steps


#### Step 1: assembly ####
assemblies = Array.new
read_1s = inputfiles  # with relative path
read_2s = []  # with relative path
sample_names = []  # without relative path
sample_name_id = {}  # without relative path
read1_pattern = ""
read2_pattern = ""
read_1s.each_with_index do |read_1, idx|
  if read_1.include?("_R1_") then
    read1_pattern = "_R1_"
    read2_pattern = "_R2_"
    read_2 = read_1.sub(read1_pattern,read2_pattern)
  elsif read_1.include?("_1") then
    read1_pattern = "_1"
    read2_pattern = "_2"
    read_2 = read_1.sub(read1_pattern,read2_pattern)
  else
    puts "The pattern of fastq.gz files is not recognized."
    exit(-1)
  end
  read_2s << read_2

  sample_name = read_1.split("/")[-1].chomp(".fastq.gz")
  sample_name_id[sample_name] = idx 
  sample_names << sample_name
  
  if steps.include?("1") then
    puts "running spades.py --careful --threads 15 -1 #{read_1} -2 #{read_2} -o #{outputdir}/#{sample_name}_assembly"
    `#{SPADES} --careful --threads 15 -1 #{read_1} -2 #{read_2} -o #{outputdir}/#{sample_name}_assembly`
    `samtools faidx #{outputdir}/#{sample_name}_assembly/contigs.fasta`
  end
  
  assembly = "#{outputdir}/#{sample_name}_assembly/contigs.fasta"
  assemblies << assembly  # with relative path
end


#### Step 2: mapping and filtering ####
# check if a readid exists in all mappings
# a read from sample j, it should all mapped on different assemblies[i]
# check if it exists on all mappings
0.upto(sample_names.length-1) do |i|
  if steps.include?("2") then
    # prepare read1 filepath and read2 filepath: reads from sample i
    read1_file = "#{outputdir}/#{sample_names[i]}_assembly/corrected/#{sample_names[i]}.fastq.00.0_0.cor.fastq.gz"
    read2_file = "#{outputdir}/#{sample_names[i]}_assembly/corrected/#{sample_names[i].sub(read1_pattern,read2_pattern)}.fastq.00.0_0.cor.fastq.gz"
    
    h = {}
    assemblies.each_with_index do |assembly, j| 
      ## Bowtie2 ##
      # TODO: using fork to process data simultaneously
      if j==0 then
	puts "#{BOWTIE2_BUILD} #{assemblies[j]} #{assemblies[j]}_index"
	`#{BOWTIE2_BUILD} #{assemblies[j]} #{assemblies[j]}_index`
	puts "#{BOWTIE2} -x #{assemblies[j]}_index -1 #{read1_file} -2 #{read2_file} --threads 15 --rg ID:#{sample_names[i]} --rg SM:#{sample_names[i]} --rg PL:illumina --rg LB:standard --very-sensitive | samtools view -bS -f0x2 - | samtools sort -n - #{outputdir}/#{sample_names[i]}-#{sample_names[j]}_sbyn"
	#`#{BOWTIE2} -x #{assemblies[j]}_index -1 #{read1_file} -2 #{read2_file} --threads 15 --rg ID:#{sample_names[i]} --rg SM:#{sample_names[i]} --rg PL:illumina --rg LB:standard --very-sensitive | samtools view -bS -f0x2 - | samtools sort -n - #{outputdir}/#{sample_names[i]}-#{sample_names[j]}_sorted`
# 	`#{BOWTIE2} -x #{assemblies[j]}_index -1 #{read1_file} -2 #{read2_file} --threads 15 --rg ID:#{sample_names[i]} --rg SM:#{sample_names[i]} --rg PL:illumina --rg LB:standard --very-sensitive | samtools view -bS - > #{outputdir}/#{sample_names[i]}-#{sample_names[j]}.bam`
	#1 `#{BOWTIE2} -x #{assemblies[j]}_index -1 #{read1_file} -2 #{read2_file} --threads 15 --rg ID:#{sample_names[i]} --rg SM:#{sample_names[i]} --rg PL:illumina --rg LB:standard --very-sensitive | samtools view -bS -f0x2 - > #{outputdir}/#{sample_names[i]}-#{sample_names[j]}.bam`
	#2 `samtools sort -n #{outputdir}/#{sample_names[i]}-#{sample_names[j]}.bam #{outputdir}/#{sample_names[i]}-#{sample_names[j]}_sbyn`
	`#{BOWTIE2} -x #{assemblies[j]}_index -1 #{read1_file} -2 #{read2_file} --threads 15 --rg ID:#{sample_names[i]} --rg SM:#{sample_names[i]} --rg PL:illumina --rg LB:standard --very-sensitive | samtools view -bS -f0x2 - | samtools sort -n - #{outputdir}/#{sample_names[i]}-#{sample_names[j]}_sbyn`
	# https://www.biostars.org/p/56246/
	# https://www.biostars.org/p/95929/
	# https://www.biostars.org/p/44681/
	# https://www.biostars.org/p/75207/
	# http://broadinstitute.github.io/picard/explain-flags.html  
	`samtools view -h #{outputdir}/#{sample_names[i]}-#{sample_names[j]}_sbyn.bam | grep -v "XS:i:" | ./foo.py | samtools view -bS - | samtools sort -n - #{outputdir}/#{sample_names[i]}-#{sample_names[j]}_sbyn2`  ## Note that the file only contains Reads-Pairs Aligned Concordantly Exactly 1 Time.
	f = IO.popen("bamToBed -i #{outputdir}/#{sample_names[i]}-#{sample_names[j]}_sbyn2.bam -bedpe | cut -f7")
	f.each_line do |line|
	  line = line.chomp
	  h[line] = 1
	end
      else
	puts "#{BOWTIE2_BUILD} #{assemblies[j]} #{assemblies[j]}_index"
	`#{BOWTIE2_BUILD} #{assemblies[j]} #{assemblies[j]}_index`
	puts "#{BOWTIE2} -x #{assemblies[j]}_index -1 #{read1_file} -2 #{read2_file} --threads 15 --rg ID:#{sample_names[i]} --rg SM:#{sample_names[i]} --rg PL:illumina --rg LB:standard --very-sensitive | samtools view -bS -f0x2 - | samtools sort -n - #{outputdir}/#{sample_names[i]}-#{sample_names[j]}_sbyn"
	#`#{BOWTIE2} -x #{assemblies[j]}_index -1 #{read1_file} -2 #{read2_file} --threads 15 --rg ID:#{sample_names[i]} --rg SM:#{sample_names[i]} --rg PL:illumina --rg LB:standard --very-sensitive | samtools view -bS -f0x2 - | samtools sort -n - #{outputdir}/#{sample_names[i]}-#{sample_names[j]}_sorted`
	`#{BOWTIE2} -x #{assemblies[j]}_index -1 #{read1_file} -2 #{read2_file} --threads 15 --rg ID:#{sample_names[i]} --rg SM:#{sample_names[i]} --rg PL:illumina --rg LB:standard --very-sensitive | samtools view -bS -f0x2 - | samtools sort -n - #{outputdir}/#{sample_names[i]}-#{sample_names[j]}_sbyn`
	`samtools view -h #{outputdir}/#{sample_names[i]}-#{sample_names[j]}_sbyn.bam | grep -v "XS:i:" | ./foo.py | samtools view -bS - | samtools sort -n - #{outputdir}/#{sample_names[i]}-#{sample_names[j]}_sbyn2`  ## Note that the file only contains Reads-Pairs Aligned Concordantly Exactly 1 Time.
	f = IO.popen("bamToBed -i #{outputdir}/#{sample_names[i]}-#{sample_names[j]}_sbyn2.bam -bedpe | cut -f7")
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
      f_filtered_sam = File.open("#{outputdir}/#{sample_names[i]}-#{sample_names[j]}_filtered.sam", "w")
      puts "IO.popen('samtools view -h #{outputdir}/#{sample_names[i]}-#{sample_names[j]}_sbyn.bam')"
      f = IO.popen("samtools view -h #{outputdir}/#{sample_names[i]}-#{sample_names[j]}_sbyn.bam")
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
      `samtools view -bS #{outputdir}/#{sample_names[i]}-#{sample_names[j]}_filtered.sam > #{outputdir}/#{sample_names[i]}-#{sample_names[j]}_filtered.bam`
      `rm #{outputdir}/#{sample_names[i]}-#{sample_names[j]}_filtered.sam`
    end
  end
end


#### Step 3: merging bam files ####
merged_bams = []
assemblies.each_with_index do |assembly, j|
  if steps.include?("3") then
    # fix 'samtools merge' ignoring @RG
    `samtools view -H #{outputdir}/#{sample_names[0]}-#{sample_names[j]}_sorted.bam | sed '/@RG/d' - > #{outputdir}/#{sample_names[j]}_header.sam`  
    filtered_bams = []
    
    0.upto(sample_names.length-1) do |i|  # read from different samples ==> marked with RG
      filtered_bams << "#{outputdir}/#{sample_names[i]}-#{sample_names[j]}_filtered.bam"
      `echo "@RG\tID:#{sample_names[i]}\tSM:#{sample_names[i]}\tPL:illumina\tLB:standard" >> #{outputdir}/#{sample_names[j]}_header.sam`
    end
    #`samtools merge -rh #{sample_names[j]}_header.sam #{sample_names[j]}_merged.bam #{filtered_bams.join(" ")}`
    if File.exist?("#{outputdir}/#{sample_names[j]}_merged.bam") then
      `rm #{outputdir}/#{sample_names[j]}_merged.bam`
    end
    `samtools merge -h #{outputdir}/#{sample_names[j]}_header.sam #{outputdir}/#{sample_names[j]}_merged.bam #{filtered_bams.join(" ")}`
    `samtools sort #{outputdir}/#{sample_names[j]}_merged.bam #{outputdir}/#{sample_names[j]}_merged_sorted`
    `samtools index #{outputdir}/#{sample_names[j]}_merged_sorted.bam`
    `rm #{outputdir}/#{sample_names[j]}_header.sam`
    `rm #{outputdir}/#{sample_names[j]}_merged.bam`
  end
  merged_bams << "#{outputdir}/#{sample_names[j]}_merged_sorted.bam"  # merged_bams are with relative path
end

## TODO of today: compare the 2 results ##
#### Step 4a: extract orthologous regions using mugsy ####
if steps.include?("4a") then
  #short_sample_names = sample_names.map {|name| name.split("-")[1]}
  assemblies_new = []
  assemblies.each_with_index do |assembly, i|  
    `cp #{assemblies[i]} #{assemblies[i].gsub("/", "_")}`
    assemblies_new << assemblies[i].gsub("/", "_")
  end
  puts assemblies_new
  #`source /home/jhuang/Tools/mugsy_x86-64-v1r2.2/mugsyenv.sh` # manual to the dir and source the file
  puts "#{MUGSY} --directory ./ --prefix mugsy_all_assembl #{assemblies_new.join(" ")} -duplications 1"
  `#{MUGSY} --directory ./ --prefix mugsy_all_assembl #{assemblies_new.join(" ")} -duplications 1`

  # filter the records in which mult!=9 and includes dup=number
  regions = {}
  temp_mafs = []
  f = File.open("mugsy_all_assembl.maf")
  is_valid = false
  temp_maf = nil
  f.each_line do |line|
    if line.start_with?("a") && line.include?("mult=#{assemblies.length}") && !line.include?("dup=") then  # line=~/mult=2/
      #puts "is_valid = true"
      is_valid = true
      temp_filename = line.match(/label=(\d*)/)[1]
      temp_maf = File.new("#{temp_filename}.maf", "w")
      
      temp_maf.puts("##maf version=1 scoring=mugsy")
      temp_mafs << temp_maf
    elsif line.start_with?("a") then
      #puts "is_valid = false"
      is_valid = false
    end
    if is_valid and line.start_with?("s") then
      temp_maf.puts(line)
      tokens = line.split
      id = tokens[1].split(".")[0]
      contig = tokens[1].split(".")[1..-1].join(".")
      if tokens[4] == "+" then
	len = tokens[3].to_i
	begin_pos = tokens[2].to_i + 1
	end_pos = begin_pos + len - 1
      end
      if tokens[4] == "-" then
	len = tokens[3].to_i
	total_len = tokens[5].to_i
	end_pos = total_len - tokens[2].to_i
	begin_pos = end_pos - len + 1
      end
      if regions.keys.include?(id) && regions[id].keys.include?(contig) then
	regions[id][contig] << Range.new(begin_pos, end_pos)
      elsif regions.keys.include?(id) then
	regions[id][contig] = []
	regions[id][contig] << Range.new(begin_pos, end_pos)
      else
	regions[id] = {}
	regions[id][contig] = []
	regions[id][contig] << Range.new(begin_pos, end_pos)
      end
    end
  end
  if !temp_maf.nil? then
    temp_maf.close
  end
#   temp_mafs.each do |temp_maf|
#     `rm #{temp_maf}`
#   end
  #pp regions

  #http://stackoverflow.com/questions/6017523/how-to-combine-overlapping-time-ranges-time-ranges-union
  def self.merge_ranges(ranges)
    ranges = ranges.sort_by {|r| r.first }
    *outages = ranges.shift
    ranges.each do |r|
      lastr = outages[-1]
      if lastr.last >= r.first - 1
	outages[-1] = lastr.first..[r.last, lastr.last].max
      else
	outages.push(r)
      end
    end
    outages
  end

  regions.each do |id, contig_ranges|
    bed_file = File.open("#{outputdir}/#{id}.bed", "w")
    contig_ranges.each do |contig, ranges|
      #puts "before merge"
      #pp ranges
      #puts "after merge"
      merged_ranges = merge_ranges(ranges)
      merged_ranges.each do |range|
	bed_file.puts("#{contig} #{range.begin} #{range.end}")    
      end
    end
    bed_file.close
    # in this case, the interval of a bed file is [x,y], e.g. 1,4 refers to the bases 1-4 
    `sort -k1,1 -k2,2n #{outputdir}/#{id}.bed > #{outputdir}/#{id}_sorted.bed`
  end
# 
#   # NODE_38_length_4871_cov_342.61_ID_75 1 4874
#   # NODE_90_length_388_cov_567.222_ID_179 1 390
# 
#   # before merge
#   # [37315..52956, 33163..46298, 4406..15514]
#   # after merge
#   # [4406..15514, 33163..52956]
  
  # cleanup
  `rm *.maf`
  `rm mugsy_all_assembl.*`
  assemblies_new.each do |assembly_new|
    `rm #{assembly_new}`
  end
end


#### Step 4b: extract orthologous regions by checking if the position covered by reads from different samples ####
if steps.include?("4b") then
  assemblies.each_with_index do |assembly, j|
    filtered_beds = []
    0.upto(sample_names.length-1) do |i|
      #bamToBed -i ERR103405_1-ERR101899_1_filtered.bam -bedpe | cut -f1
     
      # functions of bedops: bam2bed                 bam2starch_sge    convert2bed  gvf2bed     rmsk2starch  starchcat                  vcf2starch
      # bam2bed_gnuParallel     bedextract        gff2bed      gvf2starch  sam2bed      starchcluster_gnuParallel  wig2bed
      # bam2bed_sge             bedmap            gff2starch   psl2bed     sam2starch   starchcluster_sge          wig2starch
      # bam2starch              bedops            gtf2bed      psl2starch  sort-bed     unstarch
      # bam2starch_gnuParallel  closest-features  gtf2starch   rmsk2bed    starch       vcf2bed
      puts "bamToBed -i #{outputdir}/#{sample_names[i]}-#{sample_names[j]}_filtered.bam -bedpe"
      `bamToBed -i #{outputdir}/#{sample_names[i]}-#{sample_names[j]}_filtered.bam -bedpe | cut -f1,2,6 \
       | sort-bed - \
       | bedops --merge - \
       > #{outputdir}/#{sample_names[i]}-#{sample_names[j]}_filtered.bed`
      filtered_beds << "#{outputdir}/#{sample_names[i]}-#{sample_names[j]}_filtered.bed"
      #       `bamToBed -i input.bam \
#       bamToBed -i input.bam \
#     | sort-bed - \
#     | bedops --merge - \
#     > inputBam_bamToBedCollapsed.bed
    end
    puts "bedops --intersect #{filtered_beds.join(" ")} > #{outputdir}/#{sample_names[j]}.bed"
    # in this case, the interval of a bed file is (x,y], e.g. 0,4 refers to the bases 1-4 
    `sort -k1,1 -k2,2n #{outputdir}/#{sample_names[j]}.bed > #{outputdir}/#{sample_names[j]}_sorted.bed`
    filtered_beds.each do |filtered_bed|
      `rm #{filtered_bed}`
    end
  end
end


#### Step 5: SNP calling ####
merged_bams.each_with_index do |bam_filename, j|
  if steps.include?("5") then
    if !File.exist?(assemblies.at(j).sub(".fasta", ".dict")) then
      `#{CRESEQDIC} R=#{assemblies.at(j)} O=#{assemblies.at(j).sub(".fasta", ".dict")}`
    end
    #### using FREEBAYES ####
#     puts "#{FREEBAYES} --fasta-reference #{assemblies.at(j)} --targets #{sample_names.at(assembly_id)}.bed -b #{bam_filename} -p 1 -i -X -u | #{VCFFILTER} -f 'QUAL > 30' > #{rel_path}/#{sample_names.at(j)}_raw.vcf"
#     `#{FREEBAYES} --fasta-reference #{assemblies.at(j)} --targets #{sample_names.at(assembly_id)}.bed -b #{bam_filename} -p 1 -i -X -u | #{VCFFILTER} -f 'QUAL > 30' > #{rel_path}/#{sample_names.at(j)}_raw.vcf`

    ### using GATK ####
    puts "java -jar #{GATK} -glm SNP -R #{assemblies.at(j)} -T UnifiedGenotyper -I #{bam_filename} -o #{outputdir}/#{sample_names.at(j)}_raw.vcf -ploidy 1 -L #{sample_names.at(j)}.bed"
    `java -jar #{GATK} -glm SNP -R #{assemblies.at(j)} -T UnifiedGenotyper -I #{bam_filename} -o #{outputdir}/#{sample_names.at(j)}_raw.vcf -ploidy 1 -L #{sample_names.at(j)}.bed`
    puts "java -jar #{GATK} -glm SNP -R #{assemblies.at(j)} -T UnifiedGenotyper -I #{bam_filename} -o #{outputdir}/#{sample_names.at(j)}_raw.vcf -ploidy 1 -L #{sample_names.at(j)}.bed"
    `java -jar #{GATK} -glm SNP -R #{assemblies.at(j)} -T UnifiedGenotyper -I #{bam_filename} -o #{outputdir}/#{sample_names.at(j)}_raw.vcf -ploidy 1 -L #{sample_names.at(j)}.bed`
  end
end

#### Step 6: SNP filtering ####
merged_bams.each_with_index do |bam_filename, j|
  if steps.include?("6") then
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
    `java -jar #{GATK} -T VariantFiltration -R #{assemblies[j]} -o #{outputdir}/#{sample_names.at(j)}_filtered.vcf -V #{outputdir}/#{sample_names.at(j)}_raw.vcf \
    --clusterSize 3 --clusterWindowSize 10 \
    --filterExpression "QUAL < 30.0"  --filterName "LQFilter" \
    --filterExpression "QD < 10.0" --filterName "QDFilter"  \
    --filterExpression "MQ < 30.0" --filterName "MQFilter"  \
    --filterExpression "FS > 10.0" --filterName "FSFilter"  \
    --filterExpression "HaplotypeScore > 20.0" --filterName "HaplotypeScoreFilter"`
    header_line_no=`grep -n "#CHROM" #{outputdir}/#{sample_names.at(j)}_filtered.vcf | cut -d':' -f 1`
    puts "head -n #{header_line_no.chomp} #{outputdir}/#{sample_names.at(j)}_filtered.vcf > #{outputdir}/#{sample_names.at(j)}.snp_head"
    `head -n #{header_line_no.chomp} #{outputdir}/#{sample_names.at(j)}_filtered.vcf > #{outputdir}/#{sample_names.at(j)}.snp_head`
    `cat #{outputdir}/#{sample_names.at(j)}_filtered.vcf | grep PASS | cat #{outputdir}/#{sample_names.at(j)}.snp_head - > #{outputdir}/#{sample_names.at(j)}_passed.vcf`
  end
end


#### Step 7: calculating distance matrix ####
if steps.include?("7") then
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


  # generating distance matrix in the Phylip format and newick format using NJ
  #short_sample_names = sample_names.map {|name| name.split("-")[1]}
  f = File.new("final.dist", "w")
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
end