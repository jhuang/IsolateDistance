require 'optparse'
require 'optparse/time'
require 'ostruct'
require 'pp'

class IsolDistOptparser

  CODES = %w[iso-2022-jp shift_jis euc-jp utf8 binary]
  CODE_ALIASES = { "jis" => "iso-2022-jp", "sjis" => "shift_jis" }

  #
  # Return a structure describing the options.
  #
  def self.parse(args)
    # The options specified on the command line will be collected in *options*.
    # We set default values here.
    options = OpenStruct.new
    options.inputfiles = []
    #options.inputfilename = nil
    options.outputdir = ""
    options.steps = [] 
    #options.sn = false


    opts = OptionParser.new do |opts|
      opts.banner = "Usage: #{$0}.rb -i FILE1,FILE2,... -o outfile_rootname [options]\n"\
		    "Examples: #{$0} -i ./dataset_test/1-VA1290-h-rohde_S7_L001_R1_001_l500000.fastq.gz,./dataset_test/2-VA4851-h-rohde_S8_L001_R1_001_l500000.fastq.gz -o dataset_test_output -s 1\n"\
	            "Examples: #{$0} -i ./dataset1/ERR101899_1.fastq.gz,./dataset1/ERR101900_1.fastq.gz,./dataset1/ERR103394_1.fastq.gz -o dataset1_output -s 1\n"\
                    "Examples: #{$0} -i ./dataset1/ERR101899_1.fastq.gz,./dataset1/ERR101900_1.fastq.gz,\
./dataset1/ERR103394_1.fastq.gz,./dataset1/ERR103395_1.fastq.gz,./dataset1/ERR103396_1.fastq.gz,./dataset1/ERR103397_1.fastq.gz,\
./dataset1/ERR103398_1.fastq.gz,./dataset1/ERR103400_1.fastq.gz,./dataset1/ERR103401_1.fastq.gz,./dataset1/ERR103402_1.fastq.gz,\
./dataset1/ERR103403_1.fastq.gz,./dataset1/ERR103404_1.fastq.gz,./dataset1/ERR103405_1.fastq.gz,./dataset1/ERR159680_1.fastq.gz -o dataset1_output -s 5,6\n"
      opts.separator ""
      opts.separator "Specific options:"

      
      # Mandatory argument.
      opts.on("-i", "--input FILE1,FILE2,...", Array,
              "Read input from specified file") do |input|
        options.inputfiles = input
      end
#       opts.on("-i", "--input FILE",
#               "Read input from specified filename") do |inputfilename|
#         options.inputfilename = inputfilename
#       end

      opts.on("-o", "--output STRING",
              "Specify root name for output") do |outputdir|
        options.outputdir = outputdir
      end
      
      opts.on("-s", "--step 1,2,...", Array,
              "Specify the process index 1: assembly 2: mapping and filtering 3. merging bam files 4: SNP calling 5. SNP filtering 6: calculating distance matrix") do |step|
        options.steps = step
      end
         
      # Integer
#       opts.on("-k", "--kbest [HISHAPE NUMBER]",
#               "Choosing the [HISHAPE NUMBER] best hishapes for folding kinetics") do |kbest|
#           options.kbest = kbest
#       end
            
      # Boolean switch.      
#       opts.on("-s", "--sn",
#               "Restrict the folding space to strictly negative structure") do |sn|
#           options.sn = sn
#       end 
      
      opts.separator ""
      opts.separator "Common options:"

      # No argument, shows at tail.  This will print an options summary.
      # Try it and see!
      opts.on_tail("-h", "--help", "Produce help message") do
        puts opts
        exit
      end

    end

    opts.parse!(args)
    options
  end  # parse()

end  # class OptparseExample

#options = OptparseExample.parse(ARGV)
#pp options
