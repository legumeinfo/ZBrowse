#-----------------------------------------------------------
# Create annotations table for Medicago truncatula
#-----------------------------------------------------------
import csv, os, sys
# python create_annotations.py <io_directory>

# Extract subfields from a semicolon-delimited field like
# _Peptidename=Medtr1g004930.1;Name=Medtr1g004930.1;_Description=hypothetical protein;_Searchid=PAC31100295;_Pacid=31100295;ID=927273.1;_Longest=1;_Locusname=Medtr1g00493
def extract(s, tag) :
  i0 = s.find(tag)
  i1 = s.find('=', i0) + 1
  i2 = s.find(';', i1)
  if (i2 < 0) :
    return s[i1:]
  return s[i1:i2]

# Replace ASCII character codes with their character (or HTML equivalent)
# Note that the annotations file contains HTML character codes like &nbsp; and &beta;
# but we may leave these as they are.
def replace_special_characters(s) :
  s = s.replace('%20', ' ')
  s = s.replace('%21', '!')
  s = s.replace('%22', '&quot;')
  s = s.replace('%23', '#')
  s = s.replace('%24', '$')
  s = s.replace('%25', '%')
  s = s.replace('%26', '&amp;')
  s = s.replace('%27', '\'') # &apos; does not work in some browsers
  s = s.replace('%28', '(')
  s = s.replace('%29', ')')
  s = s.replace('%2A', '*')
  s = s.replace('%2B', '+')
  s = s.replace('%2C', ',')
  s = s.replace('%2D', '-')
  s = s.replace('%2E', '.')
  s = s.replace('%2F', '/')
  s = s.replace('%3A', ':')
  s = s.replace('%3B', ';')
  s = s.replace('%3C', '&lt;')
  s = s.replace('%3D', '=')
  s = s.replace('%3E', '&gt;')
  s = s.replace('%3F', '?')
  return s

# Add the description and locus name (if any) from the transcript files.
def read_phytozome_transcript(transcript_filename, descriptions) :
  csvin = csv.reader(open(transcript_filename, 'rb'), delimiter = '\t')
  for fields in csvin :
    if (fields[0].startswith('#')) :
      continue
    if (not fields[0].startswith('chr')) :
      continue

    if (fields[8].find('Description') >= 0) :
      chromosome = fields[0][3:]
      transcript_start = fields[3]
      transcript_end = fields[4]
      key = chromosome + '-' + transcript_start + '-' + transcript_end
      desc = extract(fields[8], 'Description')
      desc = replace_special_characters(desc)
      locus = extract(fields[8], 'Locusname')
      descriptions[key] = [ desc, locus ]
  return

# Create an annotations file for input to ZBrowse.
def create_annotations() :
  # for now, it assumes that all input files are in the same directory
  # and that the output file will be written there too
  dirio = sys.argv[1]
  os.chdir(dirio)

  csvin = csv.reader(open('medtr.A17_HM341.v4.0.gff3', 'rb'), delimiter = '\t')
  csvout = csv.writer(open('Medicago_truncatula_annotations.csv', 'wb'), delimiter = ',', quotechar = '"', quoting = csv.QUOTE_MINIMAL)

  my_descriptions = {}
  read_phytozome_transcript('Transcript-chr1-1..52991155.gff3', my_descriptions)
  read_phytozome_transcript('Transcript-chr2-1..45729672.gff3', my_descriptions)
  read_phytozome_transcript('Transcript-chr3-1..55515152.gff3', my_descriptions)
  read_phytozome_transcript('Transcript-chr4-1..56582383.gff3', my_descriptions)
  read_phytozome_transcript('Transcript-chr5-1..43630510.gff3', my_descriptions)
  read_phytozome_transcript('Transcript-chr6-1..35275713.gff3', my_descriptions)
  read_phytozome_transcript('Transcript-chr7-1..49172423.gff3', my_descriptions)
  read_phytozome_transcript('Transcript-chr8-1..45569985.gff3', my_descriptions)

  csvout.writerow([ 'chromosome', 'transcript_start', 'transcript_end', 'strand', 'id', 'name', 'description', 'locus' ])
  for fields in csvin :
    if (fields[0].startswith('#')) :
      continue
    if (not fields[0].startswith('chr')) :
      continue
    if (fields[2] == 'mRNA') :
      chromosome = fields[0][3:]
      transcript_start = fields[3]
      transcript_end = fields[4]
      strand = fields[6]
      id = extract(fields[8], 'ID')
      name = extract(fields[8], 'Name')
      data = [ chromosome, transcript_start, transcript_end, strand, id, name ]
      key = chromosome + '-' + transcript_start + '-' + transcript_end
      try :
        desc = my_descriptions[key]
      except :
        desc = [ '', '' ]
      csvout.writerow(data + desc)
  return

#-----------------------------------------------------------

if __name__ == '__main__':
  create_annotations()

#-----------------------------------------------------------
