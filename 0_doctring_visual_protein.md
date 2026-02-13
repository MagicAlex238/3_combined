"""
BALANCED_MARKERS DATAFRAME COLUMN REFERENCE
===========================================
PRECOMPUTED COLUMNS - NOT RECALCULATED
===========================================================

The following columns are ALREADY CALCULATED in the pipeline.
dont need to be recalculate again but directly used:

FOLD CHANGES & LOG TRANSFORMS (already computed):
    ✓ fold_change_2vs1, fold_change_3vs2, fold_change_3vs1
    ✓ log2fc_2vs1, log2fc_3vs2, log2fc_3vs1
    ✓ max_abs_log2fc
    ✓ norm_abund_contri
    → Use these directly for plotting; do not recalculate from mean_cat values

STATISTICAL SCORES (already computed):
    ✓ p_value_score (already -log10 transformed and normalized)
    ✓ significance_score
    ✓ frequency_score (already normalized)
    ✓ fc_score (already normalized from fold changes)
    → Use these directly; do not apply log transforms or normalization

AGGREGATE SCORES (already weighted & combined):
    ✓ overall_functional_score (already weighted by functional categories)
    ✓ overall_metal_score (already weighted by metal relevance)
    ✓ overall_synergy_score (already weighted by synergy detector)
    ✓ corrosion_relevance_score (already weighted for corrosion)
    ✓ niche_specific_score (already computed from pathway classification)
    ✓ database_combined_score (already combined from BRENDA data)
    ✓ combined_score (FINAL score - already combines all components)
    → Use these directly for ranking/filtering; do not recompute or re-weight

DERIVED STRUCTURES (already processed):
    ✓ functional_combi, mechanisms_combi, operational_combi (dicts already parsed)
    ✓ synergy_combi (already created from functional + metal + operational)
    ✓ inorganic_complex, organic_complex, biofilm_complex (already classified)
    → Use these directly; do not reconstruct from raw columns

RAW DATA COLUMNS (use these for new calculations if needed):
    - mean_cat1, mean_cat2, mean_cat3 (raw abundance means)
    - p_value (raw, not transformed)
    - functional_sub, functional_child (raw category assignments)
    - mechanisms_sub, mechanisms_child (raw mechanism assignments)
    - operational_sub (raw operational assignments)
    - consolidated_metals (raw metal lists)

===========================================================


IDENTIFIERS:
    Sites : str - Site location identifier
    idx : int - Unique row identifier  
    Genus : str - Taxonomic genus name
    protein_name : str - Protein name/function

ABUNDANCE MEASURES:
    mean_cat1 : float - Mean abundance in category 1
    mean_cat2 : float - Mean abundance in category 2
    mean_cat3 : float - Mean abundance in category 3
    fold_change_2vs1 : float - Fold change cat2/cat1
    log2fc_2vs1 : float - Log2 fold change cat2/cat1
    fold_change_3vs2 : float - Fold change cat3/cat2
    log2fc_3vs2 : float - Log2 fold change cat3/cat2
    fold_change_3vs1 : float - Fold change cat3/cat1
    log2fc_3vs1 : float - Log2 fold change cat3/cat1
    max_abs_log2fc : float - Maximum absolute log2 fold change
    norm_abund_contri : float - Normalized abundance contribution

STATISTICAL:
    p_value : float - Statistical significance
    significance_score : float - Pattern significance score
    real_score : float - Real score from pattern analysis

PATHWAY/ENZYME INFO:
    pathways : str - Metabolic pathway name (single string)
    reactions : str - Enzymatic reactions (can be long text)
    enzyme_class : str - EC classification
    hierarchy : str - Hierarchical classification (e.g., "Amino acid metabolism")

FUNCTIONAL CATEGORIES (hierarchy: sub → child):
    functional_sub : str - Single parent category (e.g., "organic_acid_metabolism", "iron_metabolism")
    functional_child : str - Single child term (e.g., "acetate", "ferredoxin")
    functional_combi : dict or None - Structure: {subcategory: [child_terms]}
                       Example: {'organic_acid_metabolism': ['acetate', 'succinate']}
                       Note: Cleaned to remove '<NA>' parents; may be None for rows without functional data
                       Safe to use for analysis after create_marker_groups processing

MECHANISMS (hierarchy: sub → child):
    mechanisms_sub : str - Single parent category (e.g., "biofilm_formation", "acid_production")
    mechanisms_child : str - Single child term (e.g., "attachment", "citric acid")
    mechanisms_combi : dict or None - Structure: {subcategory: [child_terms]}
                       Example: {'biofilm_formation': ['attachment', 'eps']}
                       Note: Cleaned to remove '<NA>' parents; may be None for rows without mechanism data
                       Safe to use for analysis after create_marker_groups processing

OPERATIONAL ENVIRONMENTAL FACTORS (hierarchy: sub → children):
    operational_sub : str - Single parent category (e.g., "direct_eet", "halogen_related")
    operational_combi : dict or None - Structure: {subcategory: [child_terms]}
                        Example: {'direct_eet': ['cytochrome', 'oxidoreductase', 'reductase']}
                        Contains MULTIPLE child terms per subcategory
                        Note: Cleaned to remove '<NA>' parents; may be None for rows without operational data
                        Safe to use for analysis after create_marker_groups processing

METALS:
    consolidated_metals : list[str] - Python list of chemical element/ion strings (cleaned, no None/NaN)
                          Example: ['Fe+2', 'Zn+2', 'H+', 'SO4-2']
                          Note: Elements include charge notation (e.g., 'Fe+2', 'Cu+2', 'Zn+2')
                          Normalized values include: 'As', 'Fe+2', 'Fe+3', 'H+', 'Na+', 'PO4-3', 'S', 'S2-', 
                          'S2O3-2', 'SO3-2', 'Ba', 'Ca', 'Co', 'Cu+2', 'Cd', 'Hg', 'K+', 'Mg+2', 'Ni+2', 
                          'Zn+2', 'NO3-', 'NO2-', 'SO4-2', 'Cr+3', 'Cl-', 'Pb', 'V5+', 'MoO4-2', 'Al', 
                          'F-', 'Mn+2', 'CO3-2', 'Se'
                          Handling: Safe to explode; compare with charges included

SYNERGY (from BRENDA database + computed):
    synergy_child_list : list[str] - List of BRENDA-derived child terms (cleaned, no None/NaN)
                         Example: ['acetyl', 'ferredoxin', 'succinate']
    synergy_sub_list : list[str] - List of BRENDA-derived subcategories (cleaned, no None/NaN)
                       Example: ['iron_metabolism', 'organic_acid_metabolism']
    synergy_description : str - Single descriptive string from BRENDA
                          Example: "Iron-Organic Acid Synergy (acid-enhanced Fe corrosion)"
    synergy_score : float - Score for synergy lists
    
    # Created by synergy detector (from functional/metal/operational columns):
    synergy_combi : list[str] - List combining functional_sub, operational_sub, metals
                    Example: ['organic_acid_metabolism', 'Fe+2', 'Ni+2', 'direct_eet']
                    Note: Cleaned; each element is a string
    synergy_combi_score : float - Score for synergy_combi

COMPLEXATION ANALYSIS (created by create_marker_groups):
    inorganic_complex : str - Semicolon-separated inorganic complexation terms
                        Example: "Fe+2; Zn+2; SO4-2; phosphate"
                        Includes: metals, oxides, mineral phases, oxyanions
    organic_complex : str - Semicolon-separated organic complexation terms
                      Example: "acetate; citrate; organic acid"
                      Includes: organic acids, carboxylates, fatty acids
    biofilm_complex : str - Semicolon-separated biofilm-related terms
                      Example: "biofilm; siderophore; metal chelation; eps"
                      Includes: biofilm, EPS, chelators, adhesion factors

PATHWAY CLASSIFICATION:
    pathway_classification : str - "niche-specific", "universal", or empty for mixed/unknown
    universal_pathways : str - Universal pathway name if applicable
    niche_specific_pathways : str - Niche-specific pathway name if applicable

CALCULATED SCORES (from scoring pipeline):
    overall_functional_score : float - Weighted functional category score
    overall_metal_score : float - Weighted metal relevance score
    overall_synergy_score : float - Weighted synergy score (from synergy detector)
    corrosion_relevance_score : float - Overall corrosion relevance
    niche_specific_score : float - Niche specificity score based on pathways
    database_combined_score : float - Database-derived component score
    combined_score : float - Final combined priority score (use for sorting/ranking)
    
    # Component scores:
    p_value_score : float - -log10(p_value) normalized
    frequency_score : float - Frequency component normalized
    fc_score : float - Fold change magnitude score

EXPLANATION:
    explanation : str - Human-readable explanation of why marker was prioritized

===================================================================
MARKER GROUPS (created by create_marker_groups function)
===================================================================
The function returns a dictionary of filtered DataFrames, where each DataFrame represents a specialized marker group. Use these DataFrames directly for visualization and analysis:

1.  **top_markers**: Highest-ranked markers based on `combined_score` ($\ge 75^\text{th}$ percentile), truncated to `top_count` (default 200).
2.  **significant\_markers**: Markers with `significance_score` $\ge 75^\text{th}$ percentile.
3.  **high\_[component]\_relevance**: Sub-groups based on individual score components ($\ge 75^\text{th}$ percentile), e.g., `high_metals_relevance`, `high_functional_relevance`.
4.  **mechanisms\_all**: All markers with non-empty, parsed content in `mechanisms_combi`.
5.  **consolidated\_metals**: Markers with $\ge 2$ entries in `consolidated_metals` OR with `overall_metal_score` $\ge 60^\text{th}$ percentile.
6.  **pathways\_all**: Markers with a non-empty `pathways` string AND `niche_specific_score` $\ge 50^\text{th}$ percentile.
7.  **functional\_all**: Markers with non-empty, parsed children in `functional_combi` AND `overall_functional_score` $\ge 60^\text{th}$ percentile.
8.  **inorganic/organic/biofilm\_complexes**: Markers with non-empty content in the corresponding `*_complex` column.
9.  **operational\_all**: All markers with non-empty, parsed content in `operational_combi`.
10. **synergy\_all**: Markers with $\ge 2$ valid, non-empty elements in the **list-based** `synergy_combi` AND `overall_synergy_score` $\ge 50^\text{th}$ percentile.
11. **high\_biological\_relevance**: Markers prioritized based on a composite biological score (combining functional, synergy, and niche scores) $\ge 75^\text{th}$ percentile of the non-zero scores.
12. **corrosion\_critical**: High-priority markers (top $20\%$ `combined_score`) with combined evidence: $\text{Score}_{\ge 80\%} \land ((\text{Metal}_{\text{Corr}} \land \text{Pathway}) \lor (\text{Pathway} \land \text{Synergy}) \lor (\text{Metal}_{\text{Corr}} \land \text{Synergy}))$.


The create_marker_groups() function creates specialized subsets of markers:

TOP-LEVEL GROUPS:
    - top_markers: Top N markers by combined_score (default: 200)
    - significant_markers: Markers above significance_score threshold

COMPONENT SCORE GROUPS (75th percentile by default):
    - high_metals_relevance: High overall_metal_score
    - high_functional_relevance: High overall_functional_score
    - high_synergy_relevance: High overall_synergy_score
    - high_corrosion_relevance: High corrosion_relevance_score
    - high_niche_relevance: High niche_specific_score

CATEGORY GROUPS:
    - mechanisms_all: All markers with valid mechanisms_combi dict
    - functional_all: Markers with functional_combi + high functional score (top 40%)
    - pathways_all: Markers with pathways + high niche score (top 50%)
    - operational_all: All markers with valid operational_combi dict
    - synergy_all: Markers with multi-element synergy_combi + high synergy score (top 50%)

METAL GROUPS:
    - consolidated_metals: Markers with 2+ metals OR high metal score (top 40%)

COMPLEXATION GROUPS (corrosion metals only: Cr+3, Cu+2, Fe+2/+3, Ni+2, Mn+2, Zn+2, MoO4-2):
    - inorganic_acid_complexes: Markers with inorganic complexation terms
    - organic_acid_complexes: Markers with organic acid complexation terms
    - biofilm_complexes: Markers with biofilm-related complexation terms
    
    Note: These groups are mutually non-exclusive (markers can appear in multiple)
    Classification based on child terms from all combi dicts + synergy lists

CRITICAL MARKERS:
    - high_biological_relevance: High combined functional + synergy + niche scores
    - corrosion_critical: Top 20% combined_score + (corrosion metals & pathways) OR 
                         (pathways & synergy) OR (metals & synergy)

===================================================================
DATA TYPES & STRUCTURE
===================================================================

LIST COLUMNS (safe to explode):
    - consolidated_metals : list[str]
    - synergy_child_list : list[str]
    - synergy_sub_list : list[str]
    - synergy_combi : list[str]

DICT COLUMNS (parse with caution):
    - functional_combi : dict or None - {str: list[str]}
    - mechanisms_combi : dict or None - {str: list[str]}
    - operational_combi : dict or None - {str: list[str]}
    Note: These are cleaned by create_marker_groups to remove '<NA>' parents
          May still be None for rows without data; check before accessing

SINGLE STRING COLUMNS:
    - functional_sub, functional_child
    - mechanisms_sub, mechanisms_child
    - operational_sub
    - pathways, niche_specific_pathways, universal_pathways
    - hierarchy, enzyme_class, protein_name
    - synergy_description, explanation
    - inorganic_complex, organic_complex, biofilm_complex (semicolon-separated)

NUMERIC COLUMNS (float):
    - All abundance measures (mean_cat*, fold_change*, log2fc*)
    - All score columns (overall_*, corrosion_*, combined_*, etc.)
    - p_value, significance_score, real_score

===================================================================
NULL HANDLING
===================================================================

- String columns: May contain empty strings '', 'nan', '<NA>', or None
- List columns: Empty lists [] for missing data (never None after preprocessing)
- Dict columns: None for missing data (not empty dicts)
- Numeric columns: NaN for missing values

Always check with pd.isna() or explicit None checks before processing.

===================================================================
USAGE RECOMMENDATIONS
===================================================================
BEFORE CALCULATING ANYTHING - CHECK IF IT EXISTS FIRST ⚠️

Common mistakes to avoid:
    ❌ DON'T: Calculate log2(fold_change) → Already in log2fc_* columns
    ❌ DON'T: Calculate -log10(p_value) → Already in p_value_score
    ❌ DON'T: Normalize scores 0-1 → Already normalized in *_score columns
    ❌ DON'T: Combine component scores → Already in combined_score
    ❌ DON'T: Calculate fold changes → Already in fold_change_* columns
    ❌ DON'T: Weight functional/metal scores → Already in overall_* columns

    ✓ DO: Use combined_score directly for ranking
    ✓ DO: Use log2fc_* columns directly for fold change plots
    ✓ DO: Use overall_*_score columns directly for component analysis
    ✓ DO: Check column existence before any calculation: if 'column_name' in df.columns

FOR PLOTTING:
    1. Single-string columns (functional_sub, mechanisms_sub, pathways):
       → Best for grouping, faceting, categorical plots
       
    2. List columns (consolidated_metals, synergy_child_list):
       → Explode first: df.explode('consolidated_metals')
       → Good for frequency plots, heatmaps, distributions
       
    3. Dict columns (functional_combi, mechanisms_combi, operational_combi):
       → Avoid direct plotting
       → If needed: flatten to (key, value) pairs or use parent columns instead
       
    4. Complexation columns (inorganic_complex, organic_complex, biofilm_complex):
       → Split on '; ' if counting individual terms
       → Use boolean masks for presence/absence

FOR ANALYSIS:
    - Ranking/sorting: Use combined_score (final priority) or component scores
    - Metal-specific: Explode consolidated_metals, then filter/group
    - Pathway analysis: Use pathways or niche_specific_pathways (single strings)
    - Synergy analysis: Use synergy_combi (list) or synergy_description (string)
    - Complexation: Use the three *_complex columns created by create_marker_groups

FOR MARKER GROUPS:
    - Access via: marker_groups = create_marker_groups(balanced_markers)
    - Each group is a filtered DataFrame with same columns as input
    - Groups are sorted by combined_score (descending)
    - Use group names as keys: marker_groups['top_markers']

===================================================================
CRITICAL NOTES
===================================================================

1. consolidated_metals contains IONS with charge notation (e.g., 'Fe+2', not 'Fe')
2. Combi dicts may be None; always check: if col is not None and isinstance(col, dict)
3. For metal comparisons, match exact ion forms: 'Fe+2', 'Fe+3', 'Zn+2', etc.
4. Complexation groups only include markers with corrosion-relevant metals
5. synergy_combi is a list of strings, not a dict like other combi columns
6. Always use combined_score for final ranking (not individual component scores)
"""