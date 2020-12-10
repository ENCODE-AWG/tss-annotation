table tss
"TSS regions as used by tss-annotation group"
(
string  chrom;          "Reference sequence chromosome or scaffold"
uint    chromStart;     "Start position in chromosome"
uint    chromEnd;       "End position in chromosome"
string  name;           "miRBase ID"
uint    score;          "Score scaled from expression to fit in range 0-1000"
char[1] strand;         "+ or - for strand"
uint    thickStart;     "Start of high density region"
uint    thickEnd;       "End of high density region"
uint    reserved;       "Region color"
uint    summit;         "A point of the highest expression. (Check program documentation for about how ties are handled)"
float   rawScore;       "Expression"
)
