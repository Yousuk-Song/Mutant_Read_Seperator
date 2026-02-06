#!/usr/bin/python

import pysam
import os
import sys

def print_usage():
    """사용법 및 예시 출력"""
    usage = """
    Usage:
        python mutant_seperator_full.py <input.bam> <chrom> <pos> <variant_base> <ref_base> <gene_name>

    Arguments:
        input.bam      : 분석할 원본 BAM 파일
        chrom          : 염색체 번호 (예: chr17)
        pos            : 확인하려는 genomic position (1-based)
        variant_base   : 변이 염기 (예: T)
        ref_base       : 참조 염기 (예: C)
        gene_name      : 유전자 이름 (파일명 구분용)

    Example:
        python mutant_seperator_full.py sample.bam chr17 7673776 T C TP53
    """
    print(usage)

def extract_read_info(bam_file, output_file):
    """BAM 파일에서 ReadID, CB, UMI 정보를 추출하여 TSV로 저장"""
    samfile = pysam.AlignmentFile(bam_file, "rb")
    
    with open(output_file, 'w') as f:
        f.write("ReadID\tCellBarcode\tUMI\n")
        for read in samfile:
            read_id = read.query_name
            # Cell barcode (CB tag)
            cell_barcode = read.get_tag('CB') if read.has_tag('CB') else 'NA'
            # UMI (im tag)
            if read.has_tag('im'):
                umis_raw = read.get_tag('im')
                umi = ','.join(umis_raw.split(','))
            else:
                umi = 'NA'
            
            f.write(f"{read_id}\t{cell_barcode}\t{umi}\n")
    samfile.close()

def main():
    # 인자가 부족할 경우 사용법 출력 후 종료
    if len(sys.argv) < 7:
        print_usage()
        sys.exit(1)

    # 1. 입력 인자 설정
    bam = sys.argv[1]
    chrom = sys.argv[2]
    pos = int(sys.argv[3])
    variant = sys.argv[4]
    reference = sys.argv[5]
    gene = sys.argv[6]

    # 파일명 정의
    variant_bam = bam.replace('.bam', f'.{chrom}.{pos}.{variant}.{gene}.var.bam')
    reference_bam = bam.replace('.bam', f'.{chrom}.{pos}.{reference}.{gene}.ref.bam')

    # 2. Mutant 및 Reference Read 분리 단계
    print(f"[*] Processing BAM: {bam}")
    bamfile = pysam.AlignmentFile(bam, "rb")
    variant_out = pysam.AlignmentFile(variant_bam, "wb", template=bamfile)
    reference_out = pysam.AlignmentFile(reference_bam, "wb", template=bamfile)

    # Pileup을 통해 특정 위치의 염기 확인
    for pileupcolumn in bamfile.pileup(chrom, pos - 1, pos):
        if pileupcolumn.pos == pos - 1:
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    base = pileupread.alignment.query_sequence[pileupread.query_position]
                    if base == variant:
                        variant_out.write(pileupread.alignment)
                    elif base == reference:
                        reference_out.write(pileupread.alignment)
    
    bamfile.close()
    variant_out.close()
    reference_out.close()

    # 인덱싱 (Samtools 필요)
    print("[*] Indexing new BAM files...")
    os.system(f'samtools index {variant_bam}')
    os.system(f'samtools index {reference_bam}')

    # 3. Read 정보 추출 단계 (TSV 생성)
    print("[*] Extracting CellBarcode and UMI info...")
    var_tsv = f'{variant_bam}.cell_barcodes.umi.tsv'
    ref_tsv = f'{reference_bam}.cell_barcodes.umi.tsv'
    
    extract_read_info(variant_bam, var_tsv)
    extract_read_info(reference_bam, ref_tsv)

    print(f"[#] Success! Outputs created:")
    print(f"    - {variant_bam} / {var_tsv}")
    print(f"    - {reference_bam} / {ref_tsv}")

if __name__ == "__main__":
    main()
