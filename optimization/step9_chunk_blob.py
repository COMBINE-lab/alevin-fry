"""Step 9: the work queue holds only the bytes of each meta-chunk, not a clone of the reader's whole buffer.

The reader's buffer starts at 512 KB and grows to fit the largest cell seen; it never shrinks, and every
meta-chunk pushed to the queue was a clone of the entire buffer. With four queue slots per worker and a
multi-megabyte buffer after the first big cell, that was hundreds of megabytes of mostly unused bytes
(and as many bytes copied). The meta-chunk iterator only parses the first `num_bytes`, so the clone is
now exactly that prefix. Bytes seen by the workers are unchanged.
"""
p = "/scratch1/zahra/source-code/libradicl/libradicl/src/readers.rs"
s = open(p).read()
old = "                record_context.clone(),\n                buf.clone(),\n            );\n"
new = "                record_context.clone(),\n                // only the bytes of this meta-chunk: the buffer itself keeps the capacity of the largest cell seen\n                buf[..cbytes as usize].to_vec(),\n            );\n"
n = s.count(old); assert n == 2, n
open(p, "w").write(s.replace(old, new)); print("step 9: meta-chunks carry only their own bytes (both readers)")
