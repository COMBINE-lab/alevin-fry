"""Step 10: no second copy of the reference names.

The reference names are needed: they come from the RAD header and are what maps the transcript
identifiers of the file to the transcript-to-gene map and to the names written in
quants_mat_cols.txt (the original alevin-fry does the same). The fork also cloned all of them into an
Arc for the worker threads, whose only consumer was the removed ground-truth oracle. The loaders of
the group, gene and oracle files run on the main thread and take a slice of the header's names.
"""
p = "/scratch1/zahra/source-code/alevin-fry/src/quant.rs"
s = open(p).read()
def rep(old, new, count=1):
    global s
    n = s.count(old); assert n == count, f"expected {count} of {old[:60]!r}, found {n}"
    s = s.replace(old, new)
rep("    ref_names: Arc<Vec<String>>,\n", "")
rep("    let ref_names_shared = Arc::new(hdr.ref_names.clone());\n", "    let ref_names_shared = &hdr.ref_names;\n")
rep("        let ref_names = ref_names_shared.clone();\n", "")
rep("            tid_to_gid,\n            ref_names,\n            junction_positions,\n", "            tid_to_gid,\n            junction_positions,\n")
open(p, "w").write(s); print("step 10: reference names shared by reference, worker copy removed")
