"""
BALANCED_MARKERS DATAFRAME COLUMN REFERENCE
=============================================
Groups Definitions
| Group Name                  | Short Description |
|-----------------------------|-------|
| top_markers                 | High-confidence subset of top N markers by combined score |
| significant_markers         | Statistically significant markers across risk categories |
| high_metals_relevance       | Markers with high score relevance for metal interactions |
| high_functional_relevance   | Markers with high score functional category relevance |
| high_synergy_relevance      | Markers with top score synergy relevance |
| high_pathway_relevance      | Markers most score relevant to known pathways |
| high_niche_relevance        | Markers with high score tier/niche specificity |
| high_corrosion_relevance    | High corrosion score relevance |
| mechanisms_all              | Markers associated with any corrosion mechanism |
| consolidated_metals         | Markers interacting with metallic species |
| pathways_all                | Markers with any annotated metabolic/functional pathway and requires high niche score (top 75%)|
| functional_all              | Markers assigned to any functional category (child or present) |
| inorganic_acid_complexes    | Markers predicted to form Fe-inorganic acids new column |
| organic_acid_complexes      | Markers linked to Fe-organic acid/acetates/oxalates new column |
| biofilm_indicators          | Markers linked to biofilm-formation, direct from notebook 6 section 7 abundance_biofilm only meaningful in biofilm_indicators|
| synergy_all                 | Markers with top synergy scores across dimensions |
| high_biological_relevance   | Markers with highest combined biological relevance |
| operational_all             | Markers annotated with operational/field environmental factors |
| corrosion_critical          | Markers with highest scores and critical metal/functional relevance(potential drivers) |


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

DERIVED STRUCTURES (already processed):Combi columns are partially cleaned but still contain many None children; they are not recommended for primary scoring or plotting. Prefer *_sub and *_child single-string columns.
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
    ❌ DON'T: Use column "Category", recompute each time please

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
7. Always recompute "Category" column from dict and integrate new pathways with function rename_pathway_single

Index(['Sites', 'idx', 'Genus', 'protein_name', 'mean_cat1', 'mean_cat2',
       'mean_cat3', 'descriptive_pattern', 'pattern', 'fold_change_2vs1',
       'log2fc_2vs1', 'fold_change_3vs2', 'log2fc_3vs2', 'fold_change_3vs1',
       'log2fc_3vs1', 'max_abs_log2fc', 'real_score', 'p_value',
       'significance_score', 'pathways', 'reactions', 'enzyme_class',
       'functional_sub', 'functional_child', 'synergy_child_list',
       'synergy_sub_list', 'synergy_description', 'consolidated_metals',
       'hierarchy', 'mechanisms_sub', 'mechanisms_child', 'operational_sub',
       'operational_combi', 'overall_functional_score', 'overall_metal_score',
       'overall_synergy_score', 'corrosion_relevance_score',
       'abundance_biofilm', 'pathway_classification', 'universal_pathways',
       'niche_specific_pathways', 'combined_score', 'p_value_score',
       'frequency_score', 'database_combined_score', 'niche_specific_score',
       'synergy_combi', 'synergy_combi_score', 'fc_score', 'explanation',
       'norm_abund_contri', 'Category', 'mechanisms_combi', 'functional_combi',
       'inorganic_complex', 'organic_complex'])

              Sites     idx              Genus               protein_name  mean_cat1  \
538  site_37  843582  Simplicispira      cysteine-synthase          0.0         
503  site_11  223299  Sediminibacterium  cysteine-synthase          0.0         
76   site_32  722716  Bradyrhizobium     cysteine-synthase          0.0         
14   site_9   182236  Achromobacter      thiosulfate-dehydrogenase  0.0         

     mean_cat2  mean_cat3   descriptive_pattern     pattern  fold_change_2vs1  
538  0.047404   0.000000   transition_exclusive  increasing  4741.436028       
503  0.000000   0.141542   severe_exclusive      increasing  1.000000          
76   0.000000   0.628408   severe_exclusive      increasing  1.000000          
14   0.000000   0.027468   severe_exclusive      increasing  1.000000  

log2fc_2vs1  fold_change_3vs2  log2fc_3vs2  fold_change_3vs1  \
412  0.0          12555.039335      10.0         12555.039335       
224  10.0         0.000011         -10.0         1.000000           
53   0.0          5829.301118       10.0         5829.301118        
454  0.0          5456.330560       10.0         5456.330560        

     log2fc_3vs1  max_abs_log2fc  real_score   p_value  
412  10.0         10.0            1.214082    0.000999  
224  0.0          10.0            1.028973    0.000999  
53   10.0         10.0            0.793423    0.000999  
454  10.0         10.0            0.929491    0.000999

significance_score                                         pathways  \
329  2.600152                                                              
589  2.600152            superpathway of N-acetylneuraminate degradation   
619  2.600152            aspartate superpathway                            
207  1.650152            superpathway of N-acetylneuraminate degradation   

                reactions                           enzyme_class  \
329  nan                   Acting on a sulfur group of donors.     
589  PYRUVFORMLY-RXN       Acyltransferases.                       
619  L-ASPARTATE-OXID-RXN  Acting on the CH-NH2 group of donors.   
207  PYRUVFORMLY-RXN       Acyltransferases.                       

              functional_sub functional_child  \
329  methanogenesis           coenzyme f420     
589  organic_acid_metabolism  3 oxoacyl         
619  organic_acid_metabolism  acetyl            
207  organic_acid_metabolism  3 oxoacyl         

                                                            synergy_child_list  
329  [iron sulfur cluster, surface]                                             
589  [acetyl, fermentation, formate, iron sulfur cluster, pyruvate, succinate]  
619  [acetyl, ferredoxin, fumarate, fumaric acid, succinate]                    
207  [acetyl, fermentation, formate, iron sulfur cluster, pyruvate, succinate]
synergy_sub_list  \
288  [o2_consumption, sulfur_metabolism, organic_acid_metabolism]   
247  [organic_acid_metabolism, metal_binding_chelation]             
615  [biofilm_formation, metal_binding_chelation]                   
497  [iron_metabolism, organic_acid_metabolism]                     

                                        synergy_description  \
288  Multi-pathway Synergy (3 categories)                     
247  TOC-chelation Synergy (TOC-chelate)                      
615  biofilm-chelate Synergy (biofilm-chelate-corrosion)      
497  Iron-Organic Acid Synergy (acid-enhanced Fe corrosion)   

                                                                               consolidated_metals  \
288  [Al+3, Cd+2, Cl-, Cu+, Fe+2, H+, Hg+2, K+, Na+, PO4-3, Pb+2, SO4-2, S-2, S2O3-2, SO3-2, Zn+2]   
247  [Cl-, Fe+2, H+, Mg+2, NO3-, PO4-3, SO4-2, S-2]                                                  
615  [Fe+2, H+, NO2-, SO4-2, S-2, SO3-2]                                                             
497  [Cl-, Fe+2, Fe+3, H+, Mg+2, Na+, SO4-2, S-2, S2O3-2, SO3-2]                                     

                   hierarchy     mechanisms_sub     mechanisms_child  \
288  Energy metabolism        o2_consumption     citric acid           
247  Carbohydrate metabolism  acid_production    citric acid           
615  Energy metabolism        biofilm_formation  iron sulfur cluster   
497  Energy metabolism        h2_consumption     sulphite              

    operational_sub  \
288  direct_eet       
247  direct_eet       
615  direct_eet       
497  direct_eet       
                                                                                                                                                                  operational_combi  \
288  {'direct_eet': ['oxidoreductase'], 'halogen_related': None, 'indirect_eet': None}                                                                                                                                                                      
247  {'direct_eet': ['oxidoreductase', 'reductase'], 'halogen_related': None, 'indirect_eet': None, 'microaerobic_conditions': None, 'ph_modulation': None}                                                                                                 
615  {'direct_eet': ['reductase'], 'halogen_related': None, 'indirect_eet': None, 'microaerobic_conditions': None, 'ph_modulation': None}                                                                                                                   
497  {'direct_eet': ['c type cytochrome', 'cytochrome', 'electron transfer', 'oxidase', 'oxidoreductase', 'reductase'], 'halogen_related': ['bromide'], 'indirect_eet': ['phenazine', 'quinone'], 'microaerobic_conditions': None, 'ph_modulation': None}   

     overall_functional_score  overall_metal_score  
288  7.780118                  6.5                  
247  11.040740                 4.0                  
615  9.549796                  6.0                  
497  9.463371                  6.0

     overall_synergy_score  corrosion_relevance_score  abundance_biofilm  \
452  6.0                    15.040740                 NaN                  
223  6.0                    15.040740                 NaN                  
372  5.4                    15.549796                 NaN                  
338  5.4                    15.549796                 NaN                  

    pathway_classification universal_pathways  \
452  niche-specific                             
223  niche-specific                             
372  niche-specific                             
338  niche-specific                             

                             niche_specific_pathways  combined_score  \
452  superpathway of n-acetylneuraminate degradation  1.414762         
223  superpathway of n-acetylneuraminate degradation  1.414762         
372  nan                                              1.474731         
338  nan                                              1.280762         

     p_value_score  
452  3.000434       
223  3.000434       
372  3.698970       
338  1.409369       

     frequency_score  database_combined_score  niche_specific_score  \
146  0.007874         24.457035                2.0                    
616  0.007874         24.457035                2.0                    
496  0.007874         24.368721                2.0                    
372  0.007874         24.349716                2.0                    

                             synergy_combi  synergy_combi_score  fc_score  \
146  [organic_acid_metabolism, direct_eet]  2.5                  2.0        
616  [organic_acid_metabolism, direct_eet]  2.5                  2.0        
496  [organic_acid_metabolism, direct_eet]  2.5                  2.0        
372  [methanogenesis, direct_eet]           2.0                  2.0        

                                                                                           explanation  
146  Database (24.5) | P-value (3.0) | Pattern significance (2.6) | Functions: organic_acid_metabolism  
616  Database (24.5) | P-value (3.0) | synergy_combi_score (2.5) | Functions: organic_acid_metabolism   
496  Database (24.4) | P-value (3.0) | synergy_combi_score (2.5) | Functions: organic_acid_metabolism   
372  Database (24.3) | P-value (3.7) | Pattern significance (2.8) | Functions: methanogenesis           

     norm_abund_contri  Category  \
19   0.121086          NaN         
133  0.153856          NaN         
650  0.121255          NaN         
142  0.051132          NaN         

                                                                                                                                                                                 mechanisms_combi  \
19   {'acid_production': None, 'biofilm_formation': None, 'carbon_metabolism': None, 'h2_consumption': ['sulphite'], 'iron_metabolism': None, 'metal_chelation': None}                              
133  {'acid_production': ['citric acid'], 'biofilm_formation': None, 'carbon_metabolism': None, 'h2_consumption': None, 'iron_metabolism': None, 'metal_chelation': None}                           
650  {'acid_production': None, 'biofilm_formation': ['formic acid'], 'carbon_metabolism': None, 'h2_consumption': None, 'iron_metabolism': None, 'metal_chelation': None, 'o2_consumption': None}   
142  {'acid_production': ['citric acid'], 'biofilm_formation': None, 'carbon_metabolism': None, 'h2_consumption': None, 'iron_metabolism': None, 'metal_chelation': None}                           

                                                            functional_combi  \
19   {'biofilm_formation': None, 'iron_metabolism': None, 'methanogenesis': None, 'organic_acid_metabolism': ['acetate']}     
133  {'biofilm_formation': None, 'iron_metabolism': None, 'methanogenesis': None, 'organic_acid_metabolism': ['3 oxoacyl']}   
650  {'biofilm_formation': ['alginate'], 'iron_metabolism': None, 'methanogenesis': None, 'organic_acid_metabolism': None}    
142  {'biofilm_formation': None, 'iron_metabolism': None, 'methanogenesis': None, 'organic_acid_metabolism': ['3 oxoacyl']}   

    inorganic_complex organic_complex  
19   NaN               NaN             
133  NaN               NaN             
650  NaN               NaN             
142  NaN               NaN     

"""