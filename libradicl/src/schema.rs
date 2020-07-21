/**
 * Single-cell equivalence class
 **/
 #[derive(Debug)]
 pub struct CellEQClass<'a> {
     // transcripts defining this eq. class
     pub transcripts: &'a Vec<u32>,
     // umis with multiplicities
     // the k-mer should be a k-mer class eventually
     pub umis: Vec<(u64, u32)>,
 }