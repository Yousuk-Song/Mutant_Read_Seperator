#!/usr/bin/python

import pysam
import sys

def extract_read_info(bam_file, output_file):
    # BAM 파일 열기
    samfile = pysam.AlignmentFile(bam_file, "rb")

    # 결과를 텍스트 파일로 저장
    with open(output_file, 'w') as f:
        # 파일의 첫 번째 줄에 헤더 추가
        f.write("ReadID\tCellBarcode\tUMI\n")

        # 모든 read를 순회하며 정보 추출
        for read in samfile:
            read_id = read.query_name

            # Cell barcode 추출
            cell_barcode = read.get_tag('CB') if read.has_tag('CB') else 'NA'

            # UMI 추출
            if read.has_tag('im'):
                umis_raw = read.get_tag('im')
                umis_list = umis_raw.split(',')
                umi = ','.join(umis_list)  # UMIs를 쉼표로 구분하여 문자열로 변환
            else:
                umi = 'NA'

            # 정보를 파일에 작성
            f.write(f"{read_id}\t{cell_barcode}\t{umi}\n")

    # BAM 파일 닫기
    samfile.close()

# 사용 예제

ref_bam_file = sys.argv[1]
ref_output_file = f'{ref_bam_file}.cell_barcodes.umi.tsv'
extract_read_info(ref_bam_file, ref_output_file)

var_bam_file = sys.argv[2]
var_output_file = f'{var_bam_file}.cell_barcodes.umi.tsv'
extract_read_info(var_bam_file, var_output_file)
