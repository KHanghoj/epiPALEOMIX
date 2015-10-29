from pysam.calignmentfile cimport AlignmentFile, AlignedSegment
cimport cython

@cython.boundscheck(False)
def open_bam(char *pathtobam, int countdown):
    cdef AlignmentFile samfile
    cdef AlignedSegment rec
    cdef int upperbound, lowerbound, current
    cdef int init = 1
    samfile = AlignmentFile(pathtobam, "rb")
    for rec in samfile:
        if rec.mapq < 25 or rec.is_unmapped:
            continue
        current=rec.alen
        if init:
            upperbound = current
            lowerbound = current
            init = 0
            continue
        if current > upperbound:
            upperbound = current
        elif current < lowerbound:
            lowerbound = current
        countdown -= 1
        if not countdown:
            break
    return (lowerbound, upperbound)

