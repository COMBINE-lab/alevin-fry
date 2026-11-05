"""Step 5: the chunk reader no longer spins while the work queue is full.

Both chunk readers in libradicl (`readers.rs`) waited for a free queue slot in an empty loop at
100% CPU. They now sleep 100 microseconds between checks; a chunk takes milliseconds to process,
so the workers never wait for the reader. Wall time is unchanged, CPU time stops being burnt.
"""
p = "/scratch1/zahra/source-code/libradicl/libradicl/src/readers.rs"
s = open(p).read()
old = "                while meta_chunk_queue.is_full() {}\n"
new = ("                // the workers are busy: wait a moment instead of spinning (a chunk takes milliseconds to process)\n"
       "                while meta_chunk_queue.is_full() {\n"
       "                    std::thread::sleep(std::time::Duration::from_micros(100));\n"
       "                }\n")
n = s.count(old); assert n == 2, n
open(p, "w").write(s.replace(old, new))
print("step 5: both readers sleep while the queue is full")
