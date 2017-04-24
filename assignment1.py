import mysql.connector, pysam, pybedtools

__author__ = 'Tamás Gutsohn'


class Assignment1:
    def __init__(self):
        ## Your gene of interest
        self.gene = "API5"

        self.geneinfo = ['API5', 'NM_001142931', 'chr11', 43333504, 43366080, '+', 13,
                         b'43333504,43342370,43342960,43343534,43344979,43348056,43349338,43350261,43351514,43352057,43356827,43357407,43363982,',
                         b'43333746,43342464,43343026,43343686,43345186,43348161,43349428,43350443,43351608,43352114,43356904,43357544,43366080,']

    def fetch_gene_coordinates(self, genome_reference, file_name):
        print("Connecting to UCSC to fetch data")

        ## Open connection
        cnx = mysql.connector.connect(host='genome-mysql.cse.ucsc.edu', user='genomep', passwd='password',
                                      db=genome_reference)

        ## Get cursor
        cursor = cnx.cursor()

        ## Build query fields
        query_fields = ["refGene.name2",
                        "refGene.name",
                        "refGene.chrom",
                        "refGene.txStart",
                        "refGene.txEnd",
                        "refGene.strand",
                        "refGene.exonCount",
                        "refGene.exonStarts",
                        "refGene.exonEnds"]

        ## Build query
        query = "SELECT DISTINCT %s from refGene" % ",".join(query_fields)

        ## Execute query
        cursor.execute(query)

        with open(file_name, "w") as fh:
            for row in cursor:
                fh.write(str(row) + "\n")

        ## Close cursor & connection
        cursor.close()
        cnx.close()

        print("Done fetching data")

    def get_sam_header(self, bamfile):
        samfile = pysam.AlignmentFile(bamfile, 'rb')
        print('Samfile header: ', samfile.header)
        samfile.close()

    def get_properly_paired_reads_of_gene(self, bamfile):
        samfile = pysam.AlignmentFile(bamfile, "rb")
        print('The properly paired reads of gene are: ')
        n = 0
        for read in samfile.fetch("11", self.geneinfo[3], self.geneinfo[4]):
            if read.is_proper_pair:
                print("Properly paired reads of gene:")
                n += 1
                print(read)
        print("Number of properly paired reads: {}".format(n))

    def get_gene_reads_with_indels(self, bamfile):
        bamFP = pysam.Samfile(bamfile, "rb")
        for read in bamFP:
            if (not (read.is_unmapped)):  # if it's mapped
                cigarLine = read.cigar

            for (cigarType, cigarLength) in cigarLine:
                try:
                    if(cigarType == 1):  # insertions
                        print(read)
                    elif(cigarType == 2):  # deletion
                        print(read)
                except:
                    print("Problem")

    def calculate_total_average_coverage(self, bamfile):
        bambed = pybedtools.bedtool.BedTool(bamfile)
        cov = bambed.total_coverage()
        print("Total average coverage: ", cov)

    def calculate_gene_average_coverage(self,bamfile):
        a = pybedtools.bedtool.BedTool(bamfile)
        cov = a.coverage(bg=True)
        print("Gene average coverage: ", cov)

    def get_number_mapped_reads(self, bamfile):
        flagstat = pysam.flagstat(bamfile)
        for line in flagstat.splitlines():
            if 'mapped' in line:
                print(line)
                break

    def get_gene_symbol(self):
        print("Gene symbol: ", self.geneinfo[0])

    def get_region_of_gene(self):
        print("Region of gene on {} starts at {} and ends at {}".format(self.geneinfo[2], self.geneinfo[3],
                                                                             self.geneinfo[4]))

    def get_number_of_exons(self):
        print("Number of exons: ", self.geneinfo[6])

    def print_summary(self):
        print("Print all results here")


if __name__ == '__main__':
    bamfile = 'HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam'
    print("Assignment 1")
    assignment1 = Assignment1()
    assignment1.print_summary()
    #assignment1.fetch_gene_coordinates("hg19", "API5.txt")
    assignment1.get_sam_header(bamfile)
    #assignment1.get_properly_paired_reads_of_gene(bamfile)
    #assignment1.get_gene_reads_with_indels(bamfile)
    #assignment1.calculate_total_average_coverage(bamfile)
    #assignment1.calculate_gene_average_coverage(bamfile)
    #assignment1.get_number_mapped_reads(bamfile)
    assignment1.get_gene_symbol()
    assignment1.get_region_of_gene()
    assignment1.get_number_of_exons()
