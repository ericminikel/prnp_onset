This respoitory holds code and results for our study of age of onset in genetic prion disease:

+ [source code](/src) (details below)
+ [figures & tables](/figures), including Excel files of the Supplementary Life Tables and Supplementary Duration Tables.
+ [life tables & other input data](/data)

The source code in `full_dataset_analysis.R`, together with the full, unredacted dataset that is **not** made public in this repository due to privacy concerns, is necessary to create the life tables as well as Tables S2 and S7, and Figures S2, S4 and S5. Although users cannot reproduce these parts of the analysis, the source code is made public so that you can at least see what we did there.

The source code in `public_dataset_analysis.R` accepts the life tables as input and is sufficient to reproduce Tables 1, 2, 3, S2, S3, S5, and S6, and Figures 1 and 2, S1, S3, and S6. The simulations are computationally intensive and are not run by default; to re-run these parts of the code, edit `TRUE` to `FALSE` in lines such as `regenerate_fig_2_data = FALSE`.

Tables S1 and S4 are based on literature review and are not derived computationally, thus, these are provided in the manuscript supplement but not anywhere in the source code. 

