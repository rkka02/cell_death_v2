#! python

# HDF5 compression function
import h5py
import hdf5plugin
import numpy as np
import re

def copy_and_compress_hdf5_with_attributes(input_file_path, output_file_path, compression_opts = {"compression":"gzip", "compression_opts":9}, **kwargs):
    """
    Copies the structure, data, and attributes of an HDF5 file to a new file, compressing all datasets using gzip.

    Parameters:
    - input_file_path: path to the input HDF5 file.
    - output_file_path: path where the output HDF5 file will be created.
    - compression: Type of compression to use. Default is 'gzip'.
    - compression_opts: Compression level (1-9). Higher values increase compression but also increase processing time. Default is 9.

    Note: This function does not return anything.
    """
    def copy_attributes(source, destination):
        """
        Copies attributes from the source to the destination.
        """
        for attr_name in source.attrs:
            destination.attrs[attr_name] = source.attrs[attr_name]
    
    def insert_resolution_info_for_imagej(resolution, group_out):
        """
        Add "element_size_um" attributes for each dataset
        """
        for key in group_out:
            item_out = group_out[key]
            if isinstance(item_out, h5py.Dataset):
                item_out.attrs['element_size_um'] = resolution
            elif isinstance(item_out, h5py.Group):
                insert_resolution_info_for_imagej(resolution, item_out)

    def recursively_copy_and_compress(group_in, group_out):
        """
        Recursively copies groups/datasets from the input file to the output file with compression and copies attributes.
        """
        copy_attributes(group_in, group_out)  # Copy attributes for the group
        ri_data_regex = re.compile(r'^/Data/(3D|2DMIP)/\d{6}/?$')
        for key in group_in:
            item_in = group_in[key]
            if ri_data_regex.match(item_in.name):
                # Copy dataset with compression and its attributes
                if isinstance(item_in, h5py.Dataset):
                    data = item_in[...]
                    if np.issubdtype(data.dtype, np.integer):
                        data = data.astype(np.float32)/1e4
                else:
                    # LDM-format: an experimental data format that was depracated
                    get_tile_attr =  lambda attr_name: item_in.attrs[attr_name][0]
                    is_uint8 = get_tile_attr('ScalarType')
                    if is_uint8:
                        data_type = np.uint8
                    else:
                        data_type = np.uint16
                    tile_count = get_tile_attr('NumberOfTiles')
                    data_ndim = item_in['TILE_00'].ndim
                    data_shape = list(get_tile_attr(f'DataSize{axis}') for axis in ('Z', 'Y', 'X')[3-data_ndim:])
                    data = np.zeros(data_shape, data_type)
                    for tile_idx in range(tile_count):
                        tile_path = f'TILE_{tile_idx:02d}'
                        get_tile_attr =  lambda attr_name: item_in[tile_path].attrs[attr_name][0]
                        sampling_step = get_tile_attr('SamplingStep')
                        if sampling_step != 1:
                            # what?! I don't know why... ask Tomocube
                            continue
                        offset = list(get_tile_attr(f'DataIndexOffsetPoint{axis}') for axis in ('Z', 'Y', 'X')[3-data_ndim:])
                        last_idx = list(get_tile_attr(f'DataIndexLastPoint{axis}') for axis in ('Z', 'Y', 'X')[3-data_ndim:])
                        mapping_range = tuple(slice(start,end + 1) for start, end in zip(offset, last_idx))
                        valid_data_range = tuple(slice(0,end - start + 1) for start, end in zip(offset, last_idx))
                        data[mapping_range] += item_in[tile_path][valid_data_range]
                    data = data.astype(np.float32)
                    if is_uint8:
                        min_RI = item_in.attrs['RImin'][0]
                        data /= 1e3
                        data += min_RI
                    else:
                        data /= 1e4
                dataset_out = group_out.create_dataset(key, data=data,  **compression_opts, **kwargs)
                copy_attributes(item_in, dataset_out)  # Copy attributes for the dataset
            elif isinstance(item_in, h5py.Group):
                # Create group in the output file, copy attributes, and recurse
                group_out_sub = group_out.create_group(key)
                recursively_copy_and_compress(item_in, group_out_sub)
        # Post-processing: Write resolution attribute for ImageJ viewer
        if 'ResolutionX' in group_out.attrs:
            resolution = []
            for coord in 'ZYX':
                attr_name = f'Resolution{coord}'
                if attr_name in group_out.attrs:
                    resolution.append(group_out.attrs[attr_name][0])
                else:
                    resolution.append(1.0)
            resolution = np.array(resolution)
            insert_resolution_info_for_imagej(resolution, group_out)

    with h5py.File(input_file_path, 'r') as file_in:
        with h5py.File(output_file_path, 'w') as file_out:
            recursively_copy_and_compress(file_in, file_out)
            file_out.attrs['FormatVersion'] = '1.0'

# command-line handler
import argparse
import os
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='python compress_tcf.py',
        description='Compress TCF file format that still compatible with ImageJ and napari'
    )
    parser.add_argument('filename', type = str, help='target TCF file')
    parser.add_argument('--overwrite', type = bool, default=False, help='overwrite existing TCF file')
    args = parser.parse_args()
    # set filenames
    src_fnmae = os.path.abspath(args.filename)
    dst_fname = os.path.join(os.path.dirname(src_fnmae), 'compressed_' + os.path.basename(src_fnmae))
    # generate compressed TCF file: gzip + shuffle + scaleoffset
    copy_and_compress_hdf5_with_attributes(src_fnmae, dst_fname, compression_opts={"compression":hdf5plugin.Zstd()}, shuffle = True, scaleoffset = 4)
    if args.overwrite:
        os.remove(src_fnmae)
        os.rename(dst_fname, src_fnmae)
