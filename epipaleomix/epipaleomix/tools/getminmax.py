import pysam

def main(pathtobam, countdown, readswithinrange=10000):
    init = 1
    within_count = readswithinrange
    with pysam.AlignmentFile(pathtobam, "rb") as samfile:
        for rec in samfile:
            if rec.alen< 30 or rec.mapq < 25 or rec.is_unmapped:
                continue
            current=rec.alen
            if init:
                upperbound = current
                lowerbound = current
                init = 0
                continue
            if current > upperbound:
                upperbound = current
                within_count = readswithinrange
            elif current < lowerbound:
                lowerbound = current
                within_count = readswithinrange
            else:
                within_count -= 1
            if not within_count:
                break
            countdown -= 1
            if not countdown:
                break
    return (lowerbound, upperbound)


if __name__ == '__main__':
    print main('/Users/krishang/Desktop/example/saqqaq_chrom22.bam', 500000, 50000)
