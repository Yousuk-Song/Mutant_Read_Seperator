#!/usr/bin/python

import pysam
import os
import sys


# BAM 파일 열기
bam = sys.argv[1]

bamfile = pysam.AlignmentFile(bam, "rb")
# 확인하고 싶은 염기 위치 (0-based 인덱스)
chrom = sys.argv[2]  # 염색체 이름
pos = int(sys.argv[3])    # 확인하려는 1-based 위치
variant = sys.argv[4]
reference = sys.argv[5]
gene = sys.argv[6]

variant_bam = bam.replace('.bam', f'.{chrom}.{pos}.{variant}.{gene}.var.bam')
reference_bam = bam.replace('.bam', f'.{chrom}.{pos}.{reference}.{gene}.ref.bam')

variant_bamfile = pysam.AlignmentFile(variant_bam, "wb", template = bamfile)
reference_bamfile = pysam.AlignmentFile(reference_bam, "wb", template = bamfile)

# 특정 위치의 reads에서 염기 ‘C’ 확인
for pileupcolumn in bamfile.pileup(chrom, pos - 1, pos):
    if pileupcolumn.pos == pos - 1:  # 0-based 위치이므로 -1 해줌
#        print(f"Position {pileupcolumn.pos + 1}:")
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:  # deletion과 refskip 처리
                base = pileupread.alignment.query_sequence[pileupread.query_position]
                if base == variant:  # 염기가 ‘C’인 read만 출력
#                    print(f"Read: {pileupread.alignment.query_name}, Base: {base}")
                    variant_bamfile.write(pileupread.alignment)
                elif base == reference:
                    reference_bamfile.write(pileupread.alignment)
# BAM 파일 닫기
bamfile.close()
variant_bamfile.close()
reference_bamfile.close()

os.system(f'samtools index {variant_bam}')
os.system(f'samtools index {reference_bam}')
