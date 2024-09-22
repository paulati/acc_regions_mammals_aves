line1 = "chr5_acc_50_168262496_168262546.maf	"

line2 = "chr5_acc_50_37064337_37064387.maf	A------------TGTGG-TTTTACAA----G-------GATA-TGTTAA-----------CCCCTTGGCGCTGGCGGTCGCGGTG"

tmp1 = line1.split()

tmp2 = line2.split()

print(tmp1)
print(tmp2)

if len(tmp1) == 2:
    (key, val) = tmp1
    print((key, val))
