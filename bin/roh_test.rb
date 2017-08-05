#!/usr/bin/ruby

require 'optparse'

USAGE = "Test for ROH https://en.wikipedia.org/wiki/Runs_of_Homozygosity\n\n"\
        "Usage:\n"\
        "./roh_test.rb --file=dna_raw_file.txt\n\n"

def report_a_run
  @within_run = false
  if (@genotype_in_run == '??' && @run_length > @options[:no_call_threshold]) || (@run_length > @options[:length_threshold])
    bases = (@prev_position - @run_start_position).to_f / 1000000

    chr_str = {
      23 => 'Chr X',
      24 => 'Chr Y',
      25 => 'Chr XY',
      26 => 'mtDNA'
    }[@prev_chr_number] || "Chr #{@prev_chr_number}"

    report = "#{chr_str} has a ROH of length #{@run_length} "\
          "from position #{@run_start_position} to position #{@prev_position} (%.2f Mb)" % bases

    if (@heterozygous_count > 0)
      report += "\t(#{@heterozygous_count} heterozygous SNPs treated as homozygous)"
    end

    puts report
  end
end

def homozygous?(genotype)
  genotype[0] == genotype[1]
end

@options = {
  length_threshold: 200,
  no_call_threshold: 10,
  treat_no_calls_as_homozygous: true,
  min_to_ignore_heterozygous: 150
}

# Parse params
parser = OptionParser.new do |opts|
  opts.banner = USAGE

  opts.on("-lLENGTH", "--length=LENGTH", Numeric, "Min length of ROHs to report (default: #{@options[:length_threshold]})") do |value|
    @options[:length_threshold] = value
  end

  opts.on("-nLENGTH", "--no-call-length=LENGTH", Numeric, "Min length of no-call runs to report (default: #{@options[:no_call_threshold]})") do |value|
    @options[:no_call_threshold] = value
  end

  opts.on("-t", "--[no-]treat-no-calls", TrueClass, "Treat no-calls as homozygous when finding ROHs (default: #{@options[:treat_no_calls_as_homozygous]})") do |value|
    @options[:treat_no_calls_as_homozygous] = value
  end

  opts.on("-iCOUNT", "--treat-homo=COUNT", Numeric, "Treat as homozygous any heterozygous SNP that is at least COUNT SNPs away from its nearest heterozygous SNP (default: #{@options[:max_to_ignore_heterozygous]})") do |value|
    @options[:min_to_ignore_heterozygous] = value
  end

  opts.on("-fFILE", "--file=FILE", String, "Input raw unziped DNA file") do |value|
    @options[:file] = value
  end

  opts.on("-h", "--help", "Prints this help") do
    puts opts
    exit
  end
end

begin
  parser.parse!
  raise OptionParser::MissingArgument.new("file") if @options[:file].nil?
rescue OptionParser::InvalidOption, OptionParser::MissingArgument
  $stderr.puts $!.to_s
  $stderr.puts parser
end

# Print stat
puts "File to be processed: #{@options[:file]}"
puts "ROHs of length at least #{@options[:length_threshold]} will be reported."
puts "No-call runs of length at least #{@options[:no_call_threshold]} will be reported."
puts "No-Calls will be treated as homozygous." if @options[:treat_no_calls_as_homozygous]
puts "Heterozygous SNPs that are at least #{@options[:min_to_ignore_heterozygous]} SNPs away from the nearest heterozygous SNP will be treated as homozygous."

# Detect file type
first_line = File.open(@options[:file], &:readline).downcase

@options[:file_type] = case
when first_line.include?("ancestrydna")
  'ancestrydna'
when first_line.include?("23andme")
  '23andme'
else
  'other'
end
puts "File was detected as #{@options[:file_type]}"


# Read file
File.readlines(@options[:file]).each do |line|
  next if line.strip[0] == '#'
  next if line.downcase.include?('rsid')

  if ['23andme', 'ancestrydna'].include? @options[:file_type]
    entry = line.split("\t")
    rsid = entry[0]
    chr_id = entry[1]
    position = Integer(entry[2])
    genotype = @options[:file_type] == '23andme' ? entry[3] : (entry[3] + entry[4])
  elsif line[0] == '"'
    entry = line.split(",")
    rsid = entry[0][1..entry[0].size-2]
    chr_id = entry[1][1..entry[1].size-2]
    position = entry[2][1..entry[2].size-2]
    genotype = entry[3][1..entry[3].size-2]
  else
    entry = line.split(" ")
    rsid = entry[0]
    chr_id = entry[1]
    position = entry[2]
    genotype = entry[3]
  end

  genotype = {        # genotype normalization
    '---' => '??',
    '--'  => '??',
    'CA'  => 'AC',
    'GA'  => 'AG',
    'GC'  => 'CG',
    'TA'  => 'AT',
    'TC'  => 'CT',
    'TG'  => 'GT'
  }[genotype] || genotype

  # chromosome number normalization
  chr_number = {'MT' => 26, 'X' => 23, 'Y' => 24}[chr_id] || Integer(chr_id)

  if @prev_chr_number && @prev_chr_number > chr_number
    $stderr.puts "WARNING: Chr #{chr_number} encountered after Chr #{@prev_chr_number}. The file is not properly sorted."
  end

  if @prev_chr_number == chr_number && @prev_position && @prev_position > position
    $stderr.puts "WARNING: Chr #{chr_id} position #{position} encountered after position #{@prev_position}. The file is not properly sorted."
  end

  if @prev_chr_number != chr_number && @within_run
    report_a_run
  end

  if !@options[:treat_no_calls_as_homozygous] && @within_run && genotype == '??' && @genotype_in_run != '??'
    report_a_run
  end

  if homozygous?(genotype)
    if !@within_run
      @within_run = true
      @run_start_position = position
      @run_length = 0
      @heterozygous_count = -1
      @last_heterozygous = 0
    end
    @genotype_in_run = genotype
    @run_length += 1
  elsif @within_run
    if (@run_length + 1) - @last_heterozygous > @options[:min_to_ignore_heterozygous]
      # puts "#{@run_length} - #{@last_heterozygous} = #{(@run_length + 1) - @last_heterozygous} "
      @run_length += 1
      @last_heterozygous = @run_length
      @heterozygous_count += 1
    else
      report_a_run
    end
  end

  @prev_chr_number = chr_number
  @prev_position = position
end
