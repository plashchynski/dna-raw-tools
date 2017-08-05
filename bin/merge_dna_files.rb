#!/usr/bin/ruby
#
# Usage:
# ruby ./merge_dna_files.rb file,file1,file2... > merged_file.txt
#
# Example:
# ruby ./merge_dna_files.rb AncestryDNA.txt genome_John_Doe_v4_Full_20170428065226.txt > merged_raw.txt
#
# Supports 23andMe, AncestryDNA, and Genes for Good 23andMe compatible file formats.
# The result will be in 23andMe file format.
#

def flip_val(val)
  {
    'G' => 'C',
    'T' => 'A'
  }[val] || val
end

def is_values_equal(val1, val2)
  val1, val2 = [val1, val2].map{|val| val.size == 1 ? (val + val) : val }
  val1.split('').map{|v| flip_val(v) }.sort == val2.split('').map{|v| flip_val(v) }.sort
end

def select_record(existing_record, new_record)
  if !is_values_equal(existing_record[:val], new_record[:val])
    $stderr.puts "Conflict for #{new_record[:snp_id]}: #{new_record[:file]} vs #{existing_record[:file]}, values #{new_record[:val]} vs #{existing_record[:val]}"
  end

  # Priority for 23andMe and Genes for Good formats
  selected_record = [new_record, existing_record].find {|record| ['23andMe', 'Genes for Good'].include?(record[:type]) }
  return selected_record || new_record
end
      

$snps = {}
$stat = {}
$intersections = 0

source_files = ARGV

source_files.each do |file|
  type = nil
  File.readlines(file).each do |line|
    if line.strip[0] == '#'
      type = 'ancestry' if line.include?("AncestryDNA")
      type = '23andMe' if line.include?("23andMe")
      type = 'Genes for Good' if line.include?("Genes for Good")
      next
    end

    snp = line.strip.split("\t")
    snp_id = snp[0]
    position = snp[2]
    chr = snp[1]
    val = snp[3]
    val += snp[4] if snp[4] # Ancestry values format

    next if snp_id == 'rsid' # Ancestry file format first uncommented string
    next if ['--', '00'].include?(val) # no call or void result

    # Convert Ancestry chrs to 23andMe format
    chr = {'23' => 'X', '24' => 'Y', '25' => 'X', '26' => 'MT'}[chr] || chr

    new_recod = {
      snp_id: snp_id,
      chr: chr,
      position: position,
      val: val,
      type: type,
      file: file
    }

    $stat[file] ||= 0
    $stat[file] += 1

    if existing_record = $snps[snp_id]
      $intersections += 1
      $snps[snp_id] = select_record(existing_record, new_recod)
    else
      $snps[snp_id] = new_recod
    end
  end
end

$stderr.puts "#{$intersections} intersections found."

puts "#DNA raw data file merged from:"
$stat.each do |file, snp_count|
  puts "#    #{file} with #{snp_count} SNPs"
end
puts "# #{$intersections} intersections were found among these files."
puts "#"
puts "# rsid  chromosome      position        genotype"

$snps.each do |key,values|
  puts "#{key}\t#{values[:chr]}\t#{values[:position]}\t#{values[:val]}"
end
