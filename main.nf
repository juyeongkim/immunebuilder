#!/usr/bin/env nextflow
nextflow.enable.dsl=2

FASTA = Channel.fromPath("${params.input}/*.fasta").map {[it.baseName, it]}

process ImmuneBuilder {
  tag "${sample_id}"
  publishDir "${params.output}", mode: 'copy'

  input:
  tuple val(sample_id), path(fasta_file)

  output:
  path "*.pdb"
  path "*.npy", optional: true

  script:
  """
  #!/usr/bin/env python3

  import time
  import sys
  from ImmuneBuilder import NanoBodyBuilder2
  from ImmuneBuilder import ABodyBuilder2
  from ImmuneBuilder import TCRBuilder2
  from ImmuneBuilder.util import sequence_dict_from_fasta
  import torch
  import numpy

  n_threads = ${task.cpus}
  id = "${sample_id}"

  print("Using " + str(n_threads) + " CPU threads...")
  torch.set_num_threads(n_threads)

  print("Loading fasta file...")
  seqs = sequence_dict_from_fasta(id + ".fasta")

  keys = list(seqs.keys())
  if len(keys) == 1 and keys[0] == "H":
    predictor = NanoBodyBuilder2()
  elif len(keys) == 2 and keys[0] == "H":
    predictor = ABodyBuilder2()
  elif len(keys) == 2 and keys[0] == "A":
    predictor = TCRBuilder2()
  else:
    sys.exit("Invalid FASTA file")

  print("Predicting structure...")
  res = predictor.predict(seqs)

  print("Starting refinement...")
  max_retries = ${params.max_retries}
  for attempt in range(1, max_retries + 1):
    try:
        res.save(id + ".pdb", n_threads=n_threads)
        break
    except Exception as e:
        print(f"Attempt {attempt} failed: {e}")
        if attempt == max_retries:
            print("All attempts failed. Exiting...")
            raise
        time.sleep(1)
  print("Refinement finished. Final structure saved as pdb...")

  if "${params.save_embedding}" == "true":
    embedding = torch.cat(res.encodings, dim=-1)
    numpy.save(id + ".npy", embedding.numpy())
    print("Embedding saved as npy...")
  """
}

workflow {
  ImmuneBuilder(FASTA)
}
