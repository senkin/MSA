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

def detect_and_read_csv(file_path, mutation_type, context):
    """Read a CSV file with auto-detected separator and appropriate index columns.
    
    For tab-separated files: always use index_col=0
    For comma-separated files: use index_col from get_index_columns()
    
    Args:
        file_path: Path to the file to read
        mutation_type: Type of mutation for determining index columns
        context: Context level for determining index columns
    
    Returns:
        pd.DataFrame: The loaded dataframe, or None if reading fails
    """
    try:
        # First, try to detect separator by reading first line
        with open(file_path, 'r') as f:
            first_line = f.readline().strip()
            comma_count = first_line.count(',')
            tab_count = first_line.count('\t')
        
        # Determine separator and index columns
        if tab_count > comma_count:
            print(f"Detected tab-separated file: {file_path}")
            # Tab-separated: always use index_col=0
            separator = '\t'
            index_col = 0
        else:
            print(f"Detected comma-separated file: {file_path}")
            # Comma-separated: use context-appropriate index columns
            separator = ','
            index_col = get_index_columns(mutation_type, context)
        
        # Read with determined separator and index columns
        return pd.read_csv(file_path, sep=separator, index_col=index_col)
    
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return None

def compare_index(first, second):
    """Compare indices of two dataframes and raise error if they don't match."""
    if not first.index.equals(second.index):
        print('Converted index:', first.index.to_list())
        print('Target index:', second.index.to_list())
        raise ValueError("Index mismatch, check your input data.")

def get_index_columns(mutation_type, context):
    """Determine index columns based on mutation type and context.
    
    Args:
        mutation_type: Type of mutation (e.g., 'SBS', 'DBS', 'ID')
        context: Context level (e.g., 96, 192, 288, 384, 1536, 4608)
    
    Returns:
        int or list: Index column specification for pd.read_csv
    """
    if mutation_type == 'SBS':
        if context == 96:
            return [0, 1]
        elif context in [192, 288, 384]:
            return [0, 1, 2]
        elif context in [1536, 4608]:
            return 0
        else:
            raise ValueError(f"Unsupported context {context} for SBS mutation type.")
    else:
        # Non-SBS mutation types use single index
        return 0

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

def load_signature_templates(signature_tables_path, input_signatures_prefix, mutation_types, context):
    signatures = {}
    for mutation_type in mutation_types:
        if mutation_type == 'SBS':
            index_col = get_index_columns(mutation_type, context)
            sig_file = f'{signature_tables_path}/{input_signatures_prefix}_{mutation_type}_signatures.csv'
            if context != 96:
                sig_file = f'{signature_tables_path}/{input_signatures_prefix}_{mutation_type}_{context}_signatures.csv'
            if os.path.exists(sig_file):
                signatures[mutation_type + str(context)] = pd.read_csv(sig_file, sep=',', index_col=index_col)
        else:
            index_col = 0
            sig_file = f'{signature_tables_path}/{input_signatures_prefix}_{mutation_type}_signatures.csv'
            if os.path.exists(sig_file):
                signatures[mutation_type] = pd.read_csv(sig_file, sep=',', index_col=index_col)
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

def is_already_converted(file_path, mutation_type, context, signatures):
    """Check if file is already in desired format (comma-separated with correct index).
    
    Returns:
        tuple: (is_converted, dataframe or None)
            - is_converted: True if already in desired format
            - dataframe: The loaded dataframe if successfully read, None otherwise
    """
    try:
        # Determine index columns based on mutation type and context
        index_col = get_index_columns(mutation_type, context)
        if index_col is None:
            return False, None
        
        # Try reading with auto-detected separator
        full_df = pd.read_csv(file_path, sep=None, engine='python', index_col=index_col)
        
        # Successfully read - now check if index matches expected format
        signature_key = f'{mutation_type}{context}' if mutation_type == 'SBS' else mutation_type
        
        if signature_key in signatures:
            # Check if the indices match
            if full_df.index.equals(signatures[signature_key].index):
                print(f'File {file_path} is already in desired format with correct index')
                return True, full_df
            else:
                print(f'File {file_path} has different index - will reindex')
                return False, full_df
        else:
            # No signature template to compare against - assume it needs conversion
            return False, None
            
    except Exception:
        # Read error - needs conversion
        return False, None

def process_mutation_tables(input_files_by_type, signatures, dataset_name, output_path, context, use_extractor_for_mutation_tables, specific_files_mode=False):
    """Process mutation tables (either from directory or specific files).
    
    Args:
        specific_files_mode: If True, use context parameter directly instead of inferring from filenames
    
    Returns:
        int: Number of files successfully processed/saved
    """
    make_folder_if_not_exists(f'{output_path}/{dataset_name}')
    files_processed = 0
    
    for mutation_type, files in input_files_by_type.items():
        print(f'Processing mutation type: {mutation_type}')
        print(f'Input files: {files}')
        
        for file_path in files:
            if specific_files_mode:
                # For specific files, check if already converted
                if mutation_type == 'SBS':
                    is_converted, existing_df = is_already_converted(file_path, mutation_type, context, signatures)
                else:
                    # For non-SBS, use default context value for checking
                    context_mapping = {'DBS': 78, 'ID': 83, 'SV': 32, 'CNV': 48}
                    ctx = context_mapping.get(mutation_type, 0)
                    is_converted, existing_df = is_already_converted(file_path, mutation_type, ctx, signatures)
                
                if is_converted and existing_df is not None:
                    # File is already in desired format - just save with new name
                    if mutation_type == 'SBS':
                        new_filename = f'{output_path}/{dataset_name}/WGS_{dataset_name}.{context}.csv'
                    else:
                        filename_mapping = {
                            'DBS': f'{output_path}/{dataset_name}/WGS_{dataset_name}.dinucs.csv',
                            'ID': f'{output_path}/{dataset_name}/WGS_{dataset_name}.indels.csv',
                            'SV': f'{output_path}/{dataset_name}/WGS_{dataset_name}.SV.csv',
                            'CNV': f'{output_path}/{dataset_name}/WGS_{dataset_name}.CNV.csv'
                        }
                        new_filename = filename_mapping.get(mutation_type, f'{output_path}/{dataset_name}/WGS_{dataset_name}.{mutation_type}.csv')
                    
                    existing_df.to_csv(new_filename, sep=',')
                    print(f'File already in correct format - copied to: {new_filename}')
                    files_processed += 1
                    continue
                
                # Need conversion
                # For specific files, use the mutation type and context from parameters
                if mutation_type == 'SBS':
                    mutation_type_with_context = f'SBS{context}'
                    print(f'Processing specific file: {file_path} as {mutation_type} context {context}')
                    if process_sbs_mutation_table(file_path, mutation_type_with_context, context, signatures, dataset_name, output_path):
                        files_processed += 1
                else:
                    # For non-SBS, process directly
                    context_mapping = {'DBS': '78', 'ID': '83', 'SV': '32', 'CNV': '48'}
                    mutation_type_with_context = f'{mutation_type}{context_mapping.get(mutation_type, "")}'
                    print(f'Processing specific file: {file_path} as {mutation_type}')
                    if process_non_sbs_mutation_table(file_path, mutation_type, mutation_type_with_context, signatures, dataset_name, output_path):
                        files_processed += 1
            else:
                # Original logic for directory-based inputs - infer from filename
                if use_extractor_for_mutation_tables:
                    mutation_type_with_context = Path(file_path).parent.name
                else:
                    mutation_type_with_context = Path(file_path).stem
                
                print(f'Processing file: {file_path} (context: {mutation_type_with_context})')
                
                if mutation_type == 'SBS':
                    if process_sbs_mutation_table(file_path, mutation_type_with_context, context, signatures, dataset_name, output_path):
                        files_processed += 1
                else:
                    if process_non_sbs_mutation_table(file_path, mutation_type, mutation_type_with_context, signatures, dataset_name, output_path):
                        files_processed += 1
    
    return files_processed

def process_sbs_mutation_table(file_path, mutation_type_with_context, context, signatures, dataset_name, output_path):
    """Process SBS mutation table."""
    if context == 192 and ('384' not in mutation_type_with_context or '192' not in mutation_type_with_context):
        return False
    if context != 192 and str(context) not in mutation_type_with_context:
        return False
    try:
        input_table = detect_and_read_csv(file_path, 'SBS', context)
        if input_table is None:
            return False
        
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
        return True
    except Exception as e:
        print(f'Error processing SBS mutation table {file_path}: {e}')
        return False

def process_non_sbs_mutation_table(file_path, mutation_type, mutation_type_with_context, signatures, dataset_name, output_path):
    """Process non-SBS mutation table.
    
    Returns:
        bool: True if file was successfully processed and saved, False otherwise
    """
    expected_contexts = {'DBS': '78', 'ID': '83', 'SV': '32', 'CNV': '48'}
    
    if mutation_type in expected_contexts:
        expected_context = expected_contexts[mutation_type]
        if expected_context not in mutation_type_with_context:
            return False
    try:
        print(f'Converting {mutation_type} from {file_path}')
        input_table = detect_and_read_csv(file_path, mutation_type, 0)
        if input_table is None:
            return False
        
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
        return True
    except Exception as e:
        print(f'Error processing {mutation_type} mutation table {file_path}: {e}')
        return False

def process_signature_tables(signature_files_by_type, signatures, context, output_path, reindexed_signatures_prefix, COSMIC, specific_files_mode=False):
    """Process signature tables for reindexing.

    If `specific_files_mode` is True, each provided file is checked:
      - if already comma-separated and index matches the template -> copy (rename) as-is
      - if comma-separated but index differs -> reindex (perform conversion/reindexing)
      - otherwise -> read, convert/reindex and save

    This mirrors the behavior used for mutation tables in `process_mutation_tables`.
    
    Returns:
        int: Number of files successfully processed/saved
    """
    files_processed = 0
    for mutation_type, files in signature_files_by_type.items():
        print(f'Processing signature files for mutation type: {mutation_type}')

        for signature_table_path in files:
            if specific_files_mode:
                # For specific files, check if already converted and/or reindex if needed
                if mutation_type == 'SBS':
                    # specific-file mode for SBS expects exactly one context
                    is_converted, existing_df = is_already_converted(signature_table_path, mutation_type, context, signatures)

                    if is_converted and existing_df is not None:
                        # Already comma-separated and index matches template -> save as new filename
                        filename = f'{output_path}/{reindexed_signatures_prefix}_SBS_{context}_signatures.csv'
                        if context == 96:
                            filename = filename.replace(f'_{context}', '')
                        existing_df.to_csv(filename, sep=',')
                        print(f'File already in correct format - copied to: {filename}')
                        files_processed += 1
                        continue

                    # Not already in correct format (or index differs) -> reindex via existing helper
                    if process_sbs_signature_table_with_context(signature_table_path, context, signatures, output_path, reindexed_signatures_prefix, COSMIC):
                        files_processed += 1
                else:
                    # non-SBS specific files
                    context_mapping = {'DBS': 78, 'ID': 83, 'SV': 32, 'CNV': 48}
                    expected_context = context_mapping.get(mutation_type, 0)

                    is_converted, existing_df = is_already_converted(signature_table_path, mutation_type, expected_context, signatures)

                    if is_converted and existing_df is not None:
                        filename = f'{output_path}/{reindexed_signatures_prefix}_{mutation_type}_signatures.csv'
                        existing_df.to_csv(filename, sep=',')
                        print(f'File already in correct format - copied to: {filename}')
                        files_processed += 1
                        continue

                    # Otherwise reindex / convert using existing non-SBS routine
                    if process_non_sbs_signature_table(signature_table_path, mutation_type, signatures, output_path, reindexed_signatures_prefix):
                        files_processed += 1

            else:
                # Original directory-based behavior (unchanged)
                if mutation_type == 'SBS':
                    if process_sbs_signature_table_with_context(signature_table_path, context, signatures, output_path, reindexed_signatures_prefix, COSMIC):
                        files_processed += 1
                else:
                    if process_non_sbs_signature_table(signature_table_path, mutation_type, signatures, output_path, reindexed_signatures_prefix):
                        files_processed += 1
    return files_processed

def process_sbs_signature_table_with_context(signature_table_path, context, signatures, output_path, reindexed_signatures_prefix, COSMIC):
    """Process SBS signature table with explicit context (for specific files).

    Returns:
        bool: True if processing was successful, False otherwise.
    """
    if str(context) not in signature_table_path:
        return False
    try:
        signature_table_to_reindex = detect_and_read_csv(signature_table_path, 'SBS', context)
        if signature_table_to_reindex is None:
            return False
        
        print(f'Converting signature table {signature_table_path} (SBS, context {context})')
        
        template_context = context if not COSMIC else 96
        reindexed_signatures = convert_index(signature_table_to_reindex, context=template_context)
        
        template_key = f'SBS{template_context}'
        if template_key in signatures:
            template_signatures = signatures[template_key]
            reindexed_signatures = reindexed_signatures.reindex(template_signatures.index)
            compare_index(reindexed_signatures, template_signatures)
        else:
            raise ValueError(f"No signature template found for {template_key}. Please check your signature_tables path and files.")
        
        filename = f'{output_path}/{reindexed_signatures_prefix}_SBS_{context}_signatures.csv'
        if context == 96:
            filename = filename.replace(f'_{context}', '')
        
        reindexed_signatures.to_csv(filename, sep=',')
        print(f'Saved reindexed signature table: {filename}')
        return True
    except Exception as e:
        print(f"Error processing SBS signature table {signature_table_path}: {e}")
        return False

def process_non_sbs_signature_table(signature_table_path, mutation_type, signatures, output_path, reindexed_signatures_prefix):
    """Process non-SBS signature table.
    
    Returns:
        bool: True if successfully processed, False otherwise
    """
    context_mapping = {'DBS': 78, 'ID': 83, 'SV': 32, 'CNV': 48}
    context = context_mapping.get(mutation_type)
    
    if context and str(context) not in signature_table_path:
        return False
    
    try:
        signature_table_to_reindex = detect_and_read_csv(signature_table_path, mutation_type, context)
        if signature_table_to_reindex is None:
            return False
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
        return True
    except Exception as e:
        print(f"Error processing {mutation_type} signature table {signature_table_path}: {e}")
        return False

if __name__ == '__main__':
    parser = ArgumentParser(description='Convert SigProfilerMatrixGenerator output to MSA format')
    
    # Input specification - either directory or specific file
    input_group = parser.add_mutually_exclusive_group()
    input_group.add_argument("-i", "--input_folder", dest="input_path",
                            help="Path to SigProfiler output directory (matrix generator or extractor)")
    input_group.add_argument("-I", "--input_file", dest="input_file",
                            help="Specific input file to convert")
    
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
    parser.add_argument("-c", "--context", dest="context", type=int, default=96,
                       help="SBS context (single int), e.g. 96 (default). Supported: 96, 192, 288, 1536, 4608")
    
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
    if not options.input_path and not options.input_file:
        raise ValueError("Please specify either input directory (-i) or specific input file (-I)")
    
    if not options.reindex_signatures and not options.dataset_name:
        raise ValueError("Please specify dataset name (-d) for mutation table conversion")
    
    # Validate mutation types
    supported_types = ['SBS', 'DBS', 'ID', 'SV', 'CNV']
    for mutation_type in options.mutation_types:
        if mutation_type not in supported_types:
            raise ValueError(f"Unsupported mutation type: {mutation_type}. Supported: {supported_types}")
    
    # Additional validation for specific file inputs
    if options.input_file:
        if options.reindex_signatures:
            # For signature reindexing
            if len(options.mutation_types) != 1:
                raise ValueError(f"When using specific signature file (-I with -S), please specify exactly one mutation type. Got: {options.mutation_types}")
        else:
            # For mutation table conversion
            if len(options.mutation_types) != 1:
                raise ValueError(f"When using specific mutation table file (-I), please specify exactly one mutation type. Got: {options.mutation_types}")        
        print(f"Processing specific file with mutation type: {options.mutation_types[0]}", end='')
        if options.mutation_types[0] == 'SBS':
            print(f", context: {options.context}")
        else:
            print()
    
    # Load signature templates
    signatures = load_signature_templates(
        options.signature_tables_path, 
        options.input_signatures_prefix, 
        options.mutation_types, 
        options.context
    )
    
    # Print processing information
    if options.reindex_signatures:
        signatures_type = 'COSMIC' if options.COSMIC else 'De-Novo'
        if options.input_file:
            print(f'Converting {signatures_type} signature table from specific file: {options.input_file}')
        else:
            print(f'Converting {signatures_type} signature tables from directory: {options.input_path}')
    else:
        if options.input_file:
            print(f'Converting mutation tables for dataset {options.dataset_name} from specific file: {options.input_file}')
        else:
            print(f'Converting mutation tables for dataset {options.dataset_name} from directory: {options.input_path}')
        print(f'Mutation types: {options.mutation_types}, SBS context: {options.context}')
    
    # Process files
    if options.reindex_signatures:
        if options.input_file:
            # Process specific signature files - use parameters directly
            mutation_type = options.mutation_types[0]
            signature_files_by_type = {mutation_type: [options.input_file]}
            
            # Validate that the file(s) can be read and have compatible structure
            try:
                test_table = detect_and_read_csv(options.input_file, mutation_type, options.context)
                if test_table is None:
                    raise ValueError(f"Could not read signature file {options.input_file}")
                print(f"Successfully read signature file: {options.input_file} (detected {test_table.shape[1]} signatures)")
                
                # For SBS, verify the context matches the expected number of mutation types
                if mutation_type == 'SBS':
                    context = options.context
                    expected_rows = context
                    actual_rows_full = test_table.shape[0]
                    
                    # For 192 context, input might be 384
                    if context == 192 and actual_rows_full == 384:
                        print(f"Note: Input has 384 rows, will filter out B/N to 192")
                    elif actual_rows_full != expected_rows:
                        raise ValueError(f"Warning: Expected {expected_rows} rows for context {context}, but file has {actual_rows_full} rows")
                    
            except Exception as e:
                raise ValueError(f"Error reading signature file {options.input_file}: {e}")
        else:
            # Process directory-based signature files
            signatures_type = 'COSMIC' if options.COSMIC else 'De-Novo'
            signature_files_by_type = find_signature_files_directory(
                options.input_path, options.mutation_types, signatures_type
            )
        
        files_processed = process_signature_tables(
            signature_files_by_type, signatures, options.context, 
            options.output_path, options.reindexed_signatures_prefix, options.COSMIC,
            specific_files_mode=True if options.input_file else False  # Use parameters directly, don't infer from filenames
        )
    
    else:
        # Process mutation tables
        if options.input_file:
            # Process specific mutation files - use parameters directly
            mutation_type = options.mutation_types[0]
            input_files_by_type = {mutation_type: [options.input_file]}
            
            # Validate that the file can be read and has compatible structure
            try:
                test_table = detect_and_read_csv(options.input_file, mutation_type, options.context)
                if test_table is None:
                    raise ValueError(f"Could not read mutation table file {options.input_file}")
                print(f"Successfully read mutation table: {options.input_file} (detected {test_table.shape[1]} samples)")
                
                # For SBS, verify the context
                if mutation_type == 'SBS':
                    context = options.context
                    expected_rows = context
                    actual_rows_full = test_table.shape[0]
                    
                    # For 192 context, input might be 384
                    if context == 192 and actual_rows_full == 384:
                        print(f"Note: Input has 384 rows, will filter out B/N to 192")
                    elif actual_rows_full != expected_rows:
                        raise ValueError(f"Error: Expected {expected_rows} rows for context {context}, but file has {actual_rows_full} rows")
                        
            except Exception as e:
                raise ValueError(f"Error reading mutation table {options.input_file}: {e}")
        else:
            # Process directory-based mutation files
            input_files_by_type = find_input_files_directory(
                options.input_path, options.mutation_types, options.use_extractor_for_mutation_tables
            )
        
        files_processed = process_mutation_tables(
            input_files_by_type, signatures, options.dataset_name, 
            options.output_path, options.context, options.use_extractor_for_mutation_tables,
            specific_files_mode=True if options.input_file else False  # Use parameters directly, don't infer from filenames
        )
    
    # Report completion status
    if files_processed > 0:
        print(f'Conversion completed successfully! ({files_processed} file(s) processed)')
    else:
        raise ValueError('Error: No files were successfully processed or output files created.')