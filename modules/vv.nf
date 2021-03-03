/* VV check processes
*/

process VV_RAW_READS {
  input:
    val(samples)
    path(raw_reads)
    path(vv_config)
    path(vv_output)

  output:
    path("VV_log.txt")

  script:
    """
    raw_reads_VV.py --config ${ vv_config } \
                    --samples ${ samples.join(' ') } \
                    --input ${ raw_reads } \
                    --output ${ vv_output }
    """
}
