def markers_dataframe_doc():
    """
    Reference docstring for corrosion marker DataFrames (use in notebooks as a single-source help).

    Purpose
    - Concise reference of expected columns, types, and how to handle them in analysis and plotting.
    - Targets DataFrames produced by the pipeline (pattern_data, integrated_results, classified_results,
      increasing_markers, inverse_markers, prioritized_markers, balanced_markers, group_... tables).

    High-level rules (always follow)
    - Prefer precomputed columns. If a column (e.g. p_value_score, combined_score, fc_score, frequency_score,
      database_combined_score) is present and non-null, use it. Recompute only when the precomputed column is
      entirely missing or all-null.
    - Do not overwrite upstream-computed columns unless explicitly decide and document that change.
    - For plotting: single-string categorical columns are easiest. For multi-values use lists + explode. Avoid dicts.

    CONFIRMED COLUMN LIST (name : type - short handling notes)

    IDENTIFIERS
    - Sites : str - Site location identifier. Use for grouping/faceting.
    - idx : int - Unique row identifier for traceability.
    - Genus : str - Taxonomic genus of source organism.
    - protein_name : str - Protein name or function. Use as primary label for protein-level plots.

    ABUNDANCE & EFFECT SIZE
    - mean_cat1, mean_cat2, mean_cat3 : float - Mean abundances per category. Use directly for plotting.
    - norm_abund_contri : float - Pre-normalized abundance contribution (use as-is for size/color encodings).
    - fold_change_2vs1, fold_change_3vs2, fold_change_3vs1 : float - Linear fold changes (use as-is).
    - log2fc_2vs1, log2fc_3vs2, log2fc_3vs1 : float - Log2 fold changes (use as-is).
    - max_abs_log2fc : float - Max absolute log2FC across comparisons; used to derive fc_score if fc_score missing.
    - fc_score : float - Binned effect-size score (prefer this).

    STATISTICS & SIGNIFICANCE
    - p_value : float - Raw p-value. Only convert to -log10 if p_value_score absent.
    - p_value_score : float - Precomputed -log10(p)-like score (prefer this).
    - significance_score : float - Pattern mining significance score (use as-is).
    - real_score : float - Upstream quality/confidence metric (use as-is).
    - frequency_score : float - Precomputed normalized frequency; compute from 'Frequency' only if this is missing.

    PATHWAYS & ENZYME INFO
    - pathways : str - Single pathway label. Apply rename_pathway_single before grouping/plotting if not already normalized.
    - niche_specific_pathways : str - Single pathway label for niche-specific hits; normalize similarly.
    - pathway_classification : str - "niche-specific" or "universal" (None for mixed and unknown).
    - reactions : str - Dense textual reactions; use for tables, not plots.
    - enzyme_class : str - EC class string.
    - Hierarchy : str - Hierarchical pathway classification (use for grouping/facetingpleas).

    FUNCTIONAL / MECHANISMS / OPERATIONAL (hierarchical)
    - functional_sub : str - Function subcategory (single label). Use in plots/facets.
    - functional_child : str - Specific child term (single). Use for more granular labels.
    - mechanisms_sub : str - Mechanism subcategory (single).
    - mechanisms_child : str - Mechanism child (single).
    - operational_sub : str - Operational/environmental factor (single).
    - functional_combi, mechanisms_combi, operational_combi : dict or list-of-dicts.
        - CURRENT STATUS: these "combi" columns are unreliable/corrupted in this dataset (many keys present with None children).
        - RECOMMENDATION: Do not use these combi columns for scoring or plotting until regenerated/cleaned upstream.
        - Use sub/child single-value columns instead.

    METALS (CONFIRMED)
    - consolidated_metals : list[str] - Confirmed Python list of chemical element/ion strings per row.
      Example normalized unique values (confirmed from data):
      ['As', 'Fe', 'H+', 'Na+', 'PO4-', 'S', 'S2-', 'S2O3-', 'SO3-', 'Ba', 'Ca', 'Co', 'Cu',
       'Cd', 'Hg', 'K+', 'Mg', 'Ni', 'Zn', 'NO3-', 'NO2-', 'SO4-', 'Cr', 'Cl-', 'Pb',
       'V5+', 'Mo', 'Al', 'F-', 'Mn', 'CO3-', 'Se']
      Handling:
      - This column is a true Python list (explode() is safe).
      - For plotting or metal-specific analysis: explode to one-metal-per-row, then canonicalize tokens when needed
        (metal_mapping provided globally). Do not treat consolidated_metals as a string.

    SYNERGY (BRENDA & pipeline)
    - synergy_child_list : list[str] - BRENDA-derived child terms. ploting list is ok, but after explode when needed.
    - synergy_sub_list : list[str] - BRENDA-derived subcategories. Plot after explode when needed.
    - synergy_description : str - Descriptive text from BRENDA (table-only).
    - synergy_score : float - Score for BRENDA-derived synergy 
    
    - synergy_combi : list[str] - Combined list produced by prioritize_markers (functional_sub + operational_sub + metals).
      - Use for interpretation, but canonicalize metal tokens before logic tests.
    - synergy_combi_score : float - Precomputed synergy score (prefer this).

    PATHWAY CLASSIFICATION FIELDS
    - pathway_classification : str - 'niche-specific'/'universal'/None
    - universal_pathways, niche_specific_pathways : str - single strings for display; normalize with rename_pathway_single.

    AGGREGATE & PRIORITIZATION SCORES (PREFER THESE)
    - overall_functional_score, overall_metal_score, overall_synergy_score : float - component scores.
    - corrosion_relevance_score : float - corrosion-centric weight.
    - database_combined_score : float - database-derived composite (prefer precomputed).
    - niche_specific_score : float - pathway niche emphasis (prefer precomputed).
    - combined_score : float - final prioritized ranking (use for sorting).
    - explanation : str - human-readable contributors (table-only).

    TYPES & NULLS
    - Strings: pathways, enzyme_class, hierarchy, protein_name, synergy_description, sub/child columns.
    - Floats: all score and abundance columns.
    - Lists: consolidated_metals, synergy_child_list, synergy_sub_list, synergy_combi.
    - Dict/list-of-dicts: functional_combi, mechanisms_combi, operational_combi (do not use until fixed).
    - Nulls may appear as NaN/None/empty list/None-in-dict — handle defensively.

    PARSING + HANDLING RULES (practical)
    - If 'Category' missing: create it once using Category = Sites.map(category_dict).
    - plot metals: do df_plot = df.explode('consolidated_metals').rename(columns={'consolidated_metals':'metal_code'})
      then canonicalize metal_code via metal_mapping before comparisons.
    - Always check if "combined_score" exists; if missing compute only from precomputed component columns (do not recompute components).
    - Do not use functional_combi/mechanisms_combi/operational_combi for scoring or plots until upstream combi generation is fixed.
      Temporary approach: reconstruct combi from (sub, child) pairs when child exists, but do that only if explicitly choose to.

    PLOTTING GUIDANCE (dict vs list vs single str)
    - Best practice:
      1) Single-string columns (functional_sub, mechanisms_sub, pathways) — easiest for grouping/faceting.
      2) Lists (consolidated_metals, synergy lists) — explode to one-row-per-item for frequency, heatmaps, or metal-specific violin/box plots.
      3) Dicts (functional_combi, mechanisms_combi, operational_combi) — avoid for plotting. If strictly needed, flatten to key/value pairs first.
    - Recommendation: Keep the original structured columns for provenance; create small, exploded DataFrames specifically for visualization.

    SELECTIVITY (protein vs metals) - recommended method
    - Compute using exploded metal DataFrame grouped by (protein_name, Genus):
         metal_set = set(metal_code)
         metal_count = len(metal_set)
         selectivity_class = 'selective' if metal_count <= 1 else 'promiscuous'
    - Default threshold: <= 1 considered 'selective'. Adjust threshold as needed.

    NOTE ABOUT COMBI COLUMNS (explicit)
    - No used, they are corrupted / contain many keys with None values.
    - This docstring respects that decision: treat functional_combi, mechanisms_combi, operational_combi as deprecated for analysis
      until they are regenerated/cleaned upstream.

    QUICK CHECKS 
    - Confirm consolidated_metals is list:  all(isinstance(x, list) for x in df['consolidated_metals'].fillna([]))
    - Unique metal tokens: df['consolidated_metals'].explode().unique()
    - Confirm precomputed p_value_score used: if 'p_value_score' in df.columns: use it; else compute from p_value.
    """
    pass