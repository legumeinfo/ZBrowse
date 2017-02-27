#-----------------------------------------------------------
# Create data table for Medicago truncatula
#-----------------------------------------------------------
import csv, math, os, sys
# python create_data_table.py <input_directory> <output_directory>

# Read data from a GWAS file,
# and write them to the combined GWAS data file (for input to ZBrowse).
def read_gwas_data(filename, csvout) :
  trait = filename[0:filename.find('_')]
  csvin = csv.reader(open(filename, 'rb'), delimiter = '\t')
  csvin.next() # skip header
  for fields in csvin :
    bp = fields[1]
    chromosome = fields[2][3:]
    p_value = -math.log10(float(fields[3]))
    csvout.writerow([ chromosome, bp, trait, p_value ])
  return

# Create a GWAS data file for input to ZBrowse.
def create_data_table() :
  dirout = sys.argv[2]
  os.chdir(dirout)
  csvout = csv.writer(open('Medicago_truncatula_GWAS.csv', 'wb'), delimiter = ',', quotechar = '"', quoting = csv.QUOTE_MINIMAL)
  csvout.writerow([ 'Chromosome', 'bp', 'Trait', 'P_value' ])

  dirin = sys.argv[1]
  os.chdir(dirin)
  read_gwas_data('floweringdate_results.gwas', csvout)
  read_gwas_data('height_results.gwas', csvout)
  read_gwas_data('noda_results.gwas', csvout)
  read_gwas_data('nodb_results.gwas', csvout)
  read_gwas_data('occupancyA_results.gwas', csvout)
  read_gwas_data('occupancyB_results.gwas', csvout)
  read_gwas_data('totalnod_results.gwas', csvout)
  read_gwas_data('trichomes_results.gwas', csvout)

  return

#-----------------------------------------------------------

if __name__ == '__main__':
  create_data_table()

#-----------------------------------------------------------
