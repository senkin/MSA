""" convert_SP_to_MSA.py
Module to convert mutation tables from SigProfilerMatrixGenerator or signatures
tables from SigProfilerExtractor output paths specified by -i parameter, to
comma-separated, multi-indexed tables. Output files are saved in directory
specifiable by -o flag.

Enhanced version supports specific input files via -I flag while maintaining
backward compatibility with directory-based inputs via -i flag.
"""

import os
import glob
import copy
import pandas as pd
from pathlib import Path
from argparse import ArgumentParser
from common_methods import make_folder_if_not_exists

def compare_index(first, second):
    """Compare indices of two dataframes and raise error if they don't match."""
    if not first.index.equals(second.index):
        print('Converted index:', first.index.to_list())
        print('Target index:', second.index.to_list())
        raise ValueError("Index mismatch, check your input data.")

def convert_index(input_dataframe, context=96):
    """Convert index format from SigProfiler to MSA format."""
    input_table = copy.deepcopy(input_dataframe)
    
    if context == 96:
        input_table = input_table.sort_index(level=0)
        for element in input_table.index:
            sub = element.split('[', 1)[1].split(']')[0]
            replaced_element = element.replace('['+sub+']', sub[0])
            replaced_element = sub + ':' + replaced_element
            input_table.rename(index={element: replaced_element}, inplace=True)
        input_table.index = pd.MultiIndex.from_tuples(input_table.index.str.split(':').tolist())
    elif context in [192, 384]:
        if context == 192:
            input_table = input_table[~input_table.index.str.contains("B:")]
            input_table = input_table[~input_table.index.str.contains("N:")]
        for element in input_table.index:
            sub = element.split('[', 1)[1].split(']')[0]
            replaced_element = element.replace('['+sub+']', sub[0])
            replaced_element = replaced_element.replace(':', ':' + sub + ':')
            input_table.rename(index={element: replaced_element}, inplace=True)
        input_table.index = pd.MultiIndex.from_tuples(input_table.index.str.split(':').tolist())
    elif context == 288:
        for element in input_table.index:
            sub = element.split('[', 1)[1].split(']')[0]
            replaced_element = element.replace('['+sub+']', sub[0])
            replaced_element = replaced_element.replace(':', ':' + sub + ':')
            input_table.rename(index={element: replaced_element}, inplace=True)
        input_table.index = pd.MultiIndex.from_tuples(input_table.index.str.split(':').tolist())
    
    return input_table

def load_signature_templates(signature_tables_path, input_signatures_prefix, mutation_types, contexts):
    """Load signature templates for reindexing."""
    signatures = {}
    
    for mutation_type in mutation_types:
        if mutation_type == 'SBS':
            for context in contexts:
                if context == 96:
                    index_col = [0, 1]
                    sig_file = f'{signature_tables_path}/{input_signatures_prefix}_{mutation_type}_signatures.csv'
                elif context in [192, 288]:
                    index_col = [0, 1, 2]
                    sig_file = f'{signature_tables_path}/{input_signatures_prefix}_{mutation_type}_{context}_signatures.csv'
                elif context in [1536, 4608]:
                    index_col = 0
                    sig_file = f'{signature_tables_path}/{input_signatures_prefix}_{mutation_type}_{context}_signatures.csv'
                
                if os.path.exists(sig_file):
                    signatures[mutation_type + str(context)] = pd.read_csv(sig_file, sep=',', index_col=index_col)
        else:
            sig_file = f'{signature_tables_path}/{input_signatures_prefix}_{mutation_type}_signatures.csv'
            if os.path.exists(sig_file):
                signatures[mutation_type] = pd.read_csv(sig_file, sep=',', index_col=0)
    
    return signatures

def find_input_files_directory(input_path, mutation_types, use_extractor_for_mutation_tables):
    """Find input files when using directory-based input (-i flag)."""
    input_files_by_type = {}
    
    for mutation_type in mutation_types:
        if use_extractor_for_mutation_tables:
            files = glob.glob(f'{input_path}/{mutation_type}*/Samples.txt')
        else:  # assume matrix generator output
            files = glob.glob(f'{input_path}/{mutation_type}/*{mutation_type}*')
        
        if files:
            input_files_by_type[mutation_type] = files
        else:
            print(f"Warning: No files found for mutation type {mutation_type} in {input_path}")
    
    return input_files_by_type

def find_signature_files_directory(input_path, mutation_types, signatures_type):
    """Find signature files when using directory-based input (-i flag) for signature reindexing."""
    signature_files_by_type = {}
    
    for mutation_type in mutation_types:
        files = glob.glob(f'{input_path}/{mutation_type}*/Suggested_Solution/*{signatures_type}*/Signatures/*Signatures.txt')
        if files:
            signature_files_by_type[mutation_type] = files
        else:
            print(f"Warning: No signature files found for mutation type {mutation_type} in {input_path}")
    
    return signature_files_by_type

def process_mutation_tables(input_files_by_type, signatures, dataset_name, output_path, contexts, use_extractor_for_mutation_tables):
    """Process mutation tables (either from directory or specific files)."""
    make_folder_if_not_exists(f'{output_path}/{dataset_name}')
    
    for mutation_type, files in input_files_by_type.items():
        print(f'Processing mutation type: {mutation_type}')
        print(f'Input files: {files}')
        
        for file_path in files:
            if use_extractor_for_mutation_tables:
                mutation_type_with_context = Path(file_path).parent.name
            else:
                mutation_type_with_context = Path(file_path).stem
            
            print(f'Processing file: {file_path} (context: {mutation_type_with_context})')
            
            if mutation_type == 'SBS':
                process_sbs_mutation_table(file_path, mutation_type_with_context, contexts, signatures, dataset_name, output_path)
            else:
                process_non_sbs_mutation_table(file_path, mutation_type, mutation_type_with_context, signatures, dataset_name, output_path)

def process_sbs_mutation_table(file_path, mutation_type_with_context, contexts, signatures, dataset_name, output_path):
    """Process SBS mutation table."""
    for context in contexts:
        if context == 192 and '384' not in mutation_type_with_context:
            continue
        if context != 192 and str(context) not in mutation_type_with_context:
            continue
        
        input_table = pd.read_csv(file_path, sep='\t', index_col=0)
        print(f'Converting SBS context {context} from {file_path}')
        
        input_table = convert_index(input_table, context=context)
        signature_key = f'SBS{context}'
        
        if signature_key in signatures:
            input_table = input_table.reindex(signatures[signature_key].index)
            compare_index(input_table, signatures[signature_key])
        
        new_filename = f'{output_path}/{dataset_name}/WGS_{dataset_name}.{context}.csv'
        if context == 192:
            new_filename = new_filename.replace('SBS384', 'SBS192')
        
        input_table.to_csv(new_filename, sep=',')
        print(f'Saved: {new_filename}')

def process_non_sbs_mutation_table(file_path, mutation_type, mutation_type_with_context, signatures, dataset_name, output_path):
    """Process non-SBS mutation table."""
    expected_contexts = {'DBS': '78', 'ID': '83', 'SV': '32', 'CNV': '48'}
    
    if mutation_type in expected_contexts:
        expected_context = expected_contexts[mutation_type]
        if expected_context not in mutation_type_with_context:
            return
    
    print(f'Converting {mutation_type} from {file_path}')
    input_table = pd.read_csv(file_path, sep='\t', index_col=0)
    
    if mutation_type in signatures:
        input_table.index = signatures[mutation_type].index
    
    # Generate output filename
    filename_mapping = {
        'DBS': f'{output_path}/{dataset_name}/WGS_{dataset_name}.dinucs.csv',
        'ID': f'{output_path}/{dataset_name}/WGS_{dataset_name}.indels.csv',
        'SV': f'{output_path}/{dataset_name}/WGS_{dataset_name}.SV.csv',
        'CNV': f'{output_path}/{dataset_name}/WGS_{dataset_name}.CNV.csv'
    }
    
    new_filename = filename_mapping.get(mutation_type, f'{output_path}/{dataset_name}/WGS_{dataset_name}.{mutation_type}.csv')
    input_table.to_csv(new_filename, sep=',')
    print(f'Saved: {new_filename}')

def process_signature_tables(signature_files_by_type, signatures, contexts, output_path, reindexed_signatures_prefix, COSMIC):
    """Process signature tables for reindexing."""
    for mutation_type, files in signature_files_by_type.items():
        print(f'Processing signature files for mutation type: {mutation_type}')
        
        for signature_table_path in files:
            if mutation_type == 'SBS':
                process_sbs_signature_table(signature_table_path, contexts, signatures, output_path, reindexed_signatures_prefix, COSMIC)
            else:
                process_non_sbs_signature_table(signature_table_path, mutation_type, signatures, output_path, reindexed_signatures_prefix)

def process_sbs_signature_table(signature_table_path, contexts, signatures, output_path, reindexed_signatures_prefix, COSMIC):
    """Process SBS signature table."""
    for context in contexts:
        if str(context) not in signature_table_path:
            continue
        
        signature_table_to_reindex = pd.read_csv(signature_table_path, sep='\t', index_col=0)
        print(f'Converting signature table {signature_table_path} (SBS, context {context})')
        
        template_context = context if not COSMIC else 96
        reindexed_signatures = convert_index(signature_table_to_reindex, context=template_context)
        
        template_key = f'SBS{template_context}'
        if template_key in signatures:
            template_signatures = signatures[template_key]
            reindexed_signatures = reindexed_signatures.reindex(template_signatures.index)
            compare_index(reindexed_signatures, template_signatures)
        
        filename = f'{output_path}/{reindexed_signatures_prefix}_SBS_{context}_signatures.csv'
        if context == 96:
            filename = filename.replace(f'_{context}', '')
        
        reindexed_signatures.to_csv(filename, sep=',')
        print(f'Saved reindexed signature table: {filename}')

def process_non_sbs_signature_table(signature_table_path, mutation_type, signatures, output_path, reindexed_signatures_prefix):
    """Process non-SBS signature table."""
    context_mapping = {'DBS': 78, 'ID': 83, 'SV': 32, 'CNV': 48}
    context = context_mapping.get(mutation_type)
    
    if context and str(context) not in signature_table_path:
        return
    
    signature_table_to_reindex = pd.read_csv(signature_table_path, sep='\t', index_col=0)
    print(f'Converting signature table {signature_table_path} ({mutation_type}, context {context})')
    
    # Simply overwrite index for non-SBS mutation types
    reindexed_signatures = signature_table_to_reindex
    if mutation_type in signatures:
        template_signatures = signatures[mutation_type]
        reindexed_signatures.index = template_signatures.index
        compare_index(reindexed_signatures, template_signatures)
    
    filename = f'{output_path}/{reindexed_signatures_prefix}_{mutation_type}_signatures.csv'
    reindexed_signatures.to_csv(filename, sep=',')
    print(f'Saved reindexed signature table: {filename}')

if __name__ == '__main__':
    parser = ArgumentParser(description='Convert SigProfilerMatrixGenerator output to MSA format')
    
    # Input specification - either directory or specific files
    input_group = parser.add_mutually_exclusive_group()
    input_group.add_argument("-i", "--input_folder", dest="input_path",
                            help="Path to SigProfiler output directory (matrix generator or extractor)")
    input_group.add_argument("-I", "--input_files", nargs='+', dest="input_files",
                            help="Specific input file(s) to convert")
    
    # Mode selection
    parser.add_argument("-S", "--reindex_signatures", dest="reindex_signatures", action="store_true",
                       help="Reindex signature tables (assume SP extractor output)")
    parser.add_argument("-C", "--COSMIC", dest="COSMIC", action="store_true",
                       help="Assume COSMIC signatures (forced 96 context in SBS)")
    
    # Dataset and mutation type specification
    parser.add_argument("-d", "--dataset_name", dest="dataset_name", default='',
                       help="Dataset name to use in converted filenames")
    parser.add_argument("-t", "--mutation_types", nargs='+', dest="mutation_types", default=['SBS','DBS','ID'],
                       help="Mutation types, e.g. -t SBS DBS ID (default)")
    parser.add_argument("-c", "--contexts", nargs='+', dest="contexts", type=int, default=[96, 288, 1536],
                       help="SBS contexts e.g. -c 96 288 1536 (default). Supported: 96, 192, 288, 1536, 4608")
    
    # Output and signature paths
    parser.add_argument("-o", "--output_path", dest="output_path", default='./',
                       help="Output path for converted tables (default: ./)")
    parser.add_argument("-s", "--signature_path", dest="signature_tables_path", default='signature_tables',
                       help="Path to signature tables for index templates (default: signature_tables)")
    
    # Signature prefixes
    parser.add_argument("-p", "--input_signatures_prefix", dest="input_signatures_prefix", default='sigProfiler',
                       help="Prefix in signature filenames for index templates (default: sigProfiler)")
    parser.add_argument("-n", "--reindexed_signatures_prefix", dest="reindexed_signatures_prefix", default='sigProfilerNew',
                       help="Prefix for reindexed signature filenames (default: sigProfilerNew)")
    
    # Additional options
    parser.add_argument("-E", "--use_extractor_for_mutation_tables", dest="use_extractor_for_mutation_tables", action="store_true",
                       help="Use Samples.txt from SigProfiler extractor output rather than matrix generator output")

    options = parser.parse_args()
    
    # Validate inputs
    if not options.input_path and not options.input_files:
        raise ValueError("Please specify either input directory (-i) or specific input files (-I)")
    
    if not options.reindex_signatures and not options.dataset_name:
        raise ValueError("Please specify dataset name (-d) for mutation table conversion")
    
    # Validate mutation types
    supported_types = ['SBS', 'DBS', 'ID', 'SV', 'CNV']
    for mutation_type in options.mutation_types:
        if mutation_type not in supported_types:
            raise ValueError(f"Unsupported mutation type: {mutation_type}. Supported: {supported_types}")
    
    # Load signature templates
    signatures = load_signature_templates(
        options.signature_tables_path, 
        options.input_signatures_prefix, 
        options.mutation_types, 
        options.contexts
    )
    
    # Print processing information
    if options.reindex_signatures:
        signatures_type = 'COSMIC' if options.COSMIC else 'De-Novo'
        if options.input_files:
            print(f'Converting {signatures_type} signature tables from specific files: {options.input_files}')
        else:
            print(f'Converting {signatures_type} signature tables from directory: {options.input_path}')
    else:
        if options.input_files:
            print(f'Converting mutation tables for dataset {options.dataset_name} from specific files: {options.input_files}')
        else:
            print(f'Converting mutation tables for dataset {options.dataset_name} from directory: {options.input_path}')
        print(f'Mutation types: {options.mutation_types}, SBS contexts: {options.contexts}')
    
    # Process files
    if options.reindex_signatures:
        if options.input_files:
            # Process specific signature files
            signature_files_by_type = {}
            for file_path in options.input_files:
                # Infer mutation type from filename or path
                file_stem = Path(file_path).stem.upper()
                for mut_type in options.mutation_types:
                    if mut_type in file_stem or mut_type in str(file_path).upper():
                        if mut_type not in signature_files_by_type:
                            signature_files_by_type[mut_type] = []
                        signature_files_by_type[mut_type].append(file_path)
                        break
        else:
            # Process directory-based signature files
            signatures_type = 'COSMIC' if options.COSMIC else 'De-Novo'
            signature_files_by_type = find_signature_files_directory(
                options.input_path, options.mutation_types, signatures_type
            )
        
        process_signature_tables(
            signature_files_by_type, signatures, options.contexts, 
            options.output_path, options.reindexed_signatures_prefix, options.COSMIC
        )
    
    else:
        # Process mutation tables
        if options.input_files:
            # Process specific mutation files
            input_files_by_type = {}
            for file_path in options.input_files:
                # Infer mutation type from filename or path
                file_stem = Path(file_path).stem.upper()
                for mut_type in options.mutation_types:
                    if mut_type in file_stem or mut_type in str(file_path).upper():
                        if mut_type not in input_files_by_type:
                            input_files_by_type[mut_type] = []
                        input_files_by_type[mut_type].append(file_path)
                        break
        else:
            # Process directory-based mutation files
            input_files_by_type = find_input_files_directory(
                options.input_path, options.mutation_types, options.use_extractor_for_mutation_tables
            )
        
        process_mutation_tables(
            input_files_by_type, signatures, options.dataset_name, 
            options.output_path, options.contexts, options.use_extractor_for_mutation_tables
        )
    
    print('Conversion completed successfully!')
