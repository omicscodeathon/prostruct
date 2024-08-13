from Bio.SeqUtils.ProtParam import ProteinAnalysis

# Performs basic analysis on the target sequence
# Returns the result as dictionary
# If return value is None, the operation failed
def basic_analysis(sequence):
    try:
        analysis_result = ProteinAnalysis(sequence)
        result = {
                'length': analysis_result.length,
                'amino acid count': analysis_result.count_amino_acids(),
                # 'amino acid percent': analysis_result.get_amino_acids_percent(),
                'molecular weight': analysis_result.molecular_weight(),
                'isoelectric point': analysis_result.isoelectric_point(),
                'instability index': analysis_result.instability_index(),
                'aromacity': analysis_result.aromaticity(),
                'gravy': analysis_result.gravy(),
                # 'secondary structure fraction': analysis_result.secondary_structure_fraction()
         }

        return result

    except Exception as error:
        print(error)
        return None
