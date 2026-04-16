Model with 3 classess, 5K sequences

# Data loading
DTCR.Get_Data(
    directory='...',
    aa_column_beta=...,
    count_column=...,
    v_beta_column=...,
    j_beta_column=...,
    aggregate_by_aa=True,   # collapse identical CDR3s within each sample → correct freq
)

# Evaluation: always cross-validate with 60 samples — single split is unreliable
DTCR.Monte_Carlo_CrossVal(
    test_size=0.25,          # 15 test / ~8 valid / ~37 train samples per fold
    folds=10,
    **TRAIN_KWARGS,
)



In the help description, the author of the package suggests l2_reg should be set to 0.01. Is it not too strong, comparing with your suggestions?
The author's `l2_reg=0.01` is fine for the TF version — but in the PyTorch implementation it's applied incorrectly. Let me show you why and then fix it.

In the original TF code, l2 regularization is only on the concept projection layer (the sample-level dense layer). In the PyTorch port, `weight_decay=l2_reg` is passed to `Adam` for **all** parameters — including the sequence encoder. That makes `0.01` way too strong; it crushes the convolution weights that see 300k sequences and don't need heavy regularization.

Let me find the optimizer creation in the WF Train() to fix it:

Grep "optimizer.*Adam.*weight_decay|weight_decay.*l2_reg" (in /Users/ania/Documents/learning/deeptcr/deepTCR_torch/DeepTCR.py)

7 lines of output

Read [DeepTCR.py](vscode-webview://0m11udq8qvqhclh9i2p6jd1rmjg58pr4gr922ajca26i0dh0inq5/index.html?id=5c978283-73ed-4194-8033-3f4a8bfee5f4&parentId=6&origin=63ffd1a3-8445-4331-a735-bf530318ea81&swVersion=4&extensionId=Anthropic.claude-code&platform=electron&vscode-resource-base-authority=vscode-resource.vscode-cdn.net&parentOrigin=vscode-file%3A%2F%2Fvscode-app&session=770a887a-4f5e-4050-bb9b-fd38691ee9ac#) (lines 2086-2110)
Now let me see what `_WFModel` looks like to understand which sub-modules to split:

Grep "class _WFModel" (in /Users/ania/Documents/learning/deeptcr/deepTCR_torch/DeepTCR.py)

81 lines of output

Read [DeepTCR.py](vscode-webview://0m11udq8qvqhclh9i2p6jd1rmjg58pr4gr922ajca26i0dh0inq5/index.html?id=5c978283-73ed-4194-8033-3f4a8bfee5f4&parentId=6&origin=63ffd1a3-8445-4331-a735-bf530318ea81&swVersion=4&extensionId=Anthropic.claude-code&platform=electron&vscode-resource-base-authority=vscode-resource.vscode-cdn.net&parentOrigin=vscode-file%3A%2F%2Fvscode-app&session=770a887a-4f5e-4050-bb9b-fd38691ee9ac#) (lines 348-377)

Now I have a clear picture. The fix: use two parameter groups in the optimizer — `conv_model` (encoder) gets `weight_decay=0`, while `concept_proj`, `agg_fc`, and `head` get `weight_decay=l2_reg`. Also do the same for `_SSModel` (which has the same bug). Let me check the SS optimizer:

Grep "optimizer.*Adam.*weight_decay" (in /Users/ania/Documents/learning/deeptcr/deepTCR_torch/DeepTCR.py)

1 line of output

Grep "weight_decay" (in /Users/ania/Documents/learning/deeptcr/deepTCR_torch/DeepTCR.py)

7 lines of output

Only one occurrence — WF only. SS doesn't use l2_reg at all. Now let me make the fix:
Now let me also handle the edge case where `quantitative_agg=True` adds `c_b` parameter — that's a sample-level parameter too, so it should be in `sample_params`. Let me check if `c_b` would be captured by the current groups:

The `c_b` is an `nn.Parameter` directly on the model (not inside `concept_proj`, `agg_fc`, or `head`), so it'd fall through. Let me fix that:
This is cleaner and future-proof — encoder gets no decay, everything else (concept_proj, agg_fc, head, c_b) gets `l2_reg`.

So to answer your question directly:

**`l2_reg=0.01` is correct for the TF version, but was too strong in the PyTorch port** — because PyTorch was applying it to all ~45k parameters including the convolutional encoder, while TF only applies it to the ~800-parameter concept projection layer. With the fix above, both implementations now behave the same way, and `l2_reg=0.01` is appropriate for both.

For your 3-class, 60-sample, ~300k sequence dataset, the original recommendation stands:

- Start with `l2_reg=0.01` (now correctly scoped to the sample-level layers only)
- The conv encoder is unregularized, which is appropriate since it's trained across the full 300k sequences and doesn't overfit the sample-level labels directly