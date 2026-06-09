# Semantic Retrieval Settings

This file documents the semantic retrieval configuration used by MedNLP-Hub for variable-level vector matching in the frozen reproducibility package.

## Implementation Location

The semantic retrieval implementation is in `app.py`.

Key functions and settings:

- Local embedding cache: `.embedding_cache/variable_embeddings.sqlite3` and `.embedding_cache/variable_embeddings_*.npz`.
- Local semantic model cache: a user-specified local embedding model directory.
- Embedding loading and cache lookup: `_cache_conn` and `_load_or_build_embeddings`.
- Corpus embedding precomputation: `_precompute_embeddings`.
- Vector search and ranking: `_vector_match_records`.
- Filtering entry point: `_filtered_records`.

## Data Field Used for Semantic Matching

Semantic matching is applied to the normalized `variable_list` field.

In the source workbook, this field may come from:

- `Variable_list`
- `VariableList`
- `VariableList_mapped`

During record loading, workbook columns are normalized using `COLUMN_MAP` in `app.py`. Records without a non-empty `variable_list` are excluded from the vector-search embedding index.

## Embedding Model Configuration

The frozen reproducibility package is documented as a local semantic-model workflow. The semantic model is read from a local embedding model directory rather than requiring a remote embedding service during independent reproduction.

The exact local model can be selected by the user or deployment environment. The archive should therefore record the local model identifier and model directory used for the frozen run, but the retrieval logic is model-agnostic as long as the model converts query text and record text into comparable dense vectors.

```text
SEMANTIC_MODEL_SOURCE=local
SEMANTIC_MODEL_ID=<local embedding model identifier>
SEMANTIC_MODEL_DIR=<path to local model directory included in the archive or deployment>
EMBEDDING_CACHE_SQLITE=.embedding_cache/variable_embeddings.sqlite3
EMBEDDING_CACHE_NPZ=.embedding_cache/variable_embeddings_*.npz
VECTOR_TOP_K=50
```

No private API key is needed for reproducing the frozen semantic retrieval results when the local model files and archived embedding cache are included. If a different local embedding model is selected, the record embeddings should be regenerated from the same frozen `variable_list` values, and the model identifier should be recorded in the archive metadata.

## Query Processing

For advanced search, the variable query is passed through the `variable_q` parameter.

Default behavior:

- `variable_match_mode=vector`
- `limit=50` for variable-level semantic matching in `_filtered_records`
- results sorted by `vector_score` when no other sort mode is selected

The query text is stripped and lowercased for fallback lexical matching. For vector matching, the original stripped query string is encoded by the local semantic model or matched against the archived embedding cache.

## Embedding Normalization and Similarity

All embedding vectors are converted to `float32` arrays.

Each vector is L2-normalized:

```text
normalized_vector = vector / L2_norm(vector)
```

Similarity is computed as the dot product between the normalized query embedding and each normalized record embedding:

```text
vector_score = query_embedding dot record_embedding
```

Because both vectors are normalized, this score is cosine similarity.

## Ranking and Cutoff

Records are ranked by descending cosine similarity.

The default maximum number of returned vector-matched records is:

```text
VECTOR_TOP_K=50
```

In the current advanced-search route, variable-level vector matching explicitly uses:

```text
limit=50
```

No fixed cosine-similarity threshold is applied in the application code. The top-ranked records are returned up to the configured limit.

## Local Cache and Reproducibility

Record-level embeddings are cached in:

```text
.embedding_cache/variable_embeddings.sqlite3
.embedding_cache/variable_embeddings_*.npz
```

The cache key uses:

- the embedding model name
- SHA-256 hash of the source text
- the original text value

Cached vectors are reused only when the stored text exactly matches the current `variable_list` value and the vector dimension matches the stored metadata. This prevents stale embeddings from being used after text edits.

For full reproducibility, the frozen archive should include:

- the local embedding model files or a clear identifier/version for the local model used;
- the cached SQLite embedding database used for the benchmark export;
- the cached `.npz` vector file if it was used during export;
- the frozen `database_FINAL.xlsx` file from which the `variable_list` values were read;
- the full returned candidate lists produced from the frozen local setup.

## Fallback Behavior

If local embeddings are unavailable, MedNLP-Hub can fall back to case-insensitive substring matching against `variable_list`.

This fallback keeps the web interface usable but should not be used for the reported semantic-retrieval benchmark unless explicitly stated. The benchmark archive should therefore include a log or metadata note confirming that vector embeddings were available during candidate-list export.
