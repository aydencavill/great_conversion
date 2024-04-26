import shutil
import sys
import json
import re
import json
import zarr
import subprocess
import argparse
import tempfile
import numpy as np
from collections import defaultdict
from PyReconstruct.modules.datatypes import Series
from PyReconstruct.modules.backend.view.trace_layer import hashName


temp_jser = tempfile.NamedTemporaryFile(suffix=".jser").name
print(temp_jser)

parser = argparse.ArgumentParser()

parser.add_argument("input_jser", type=str, nargs="?", help="Filepath of a valid jser file.")

parser.add_argument(
    "--output",
    "-o",
    type=str,
    default='pyrecon_series.zarr',
    help="Optional output path",
)

parser.add_argument(
    "--groups",
    "-g",
    type=str,
    action="append",
    nargs="*",
    default=None)

args = parser.parse_args()

shutil.copy(args.input_jser,temp_jser)

run_command = ['ng-create-zarr',temp_jser]
if args.output:
    run_command += ['-o',args.output]
if args.groups:
    print(args.groups)
    groups = " ".join(args.groups[0])
    run_command += ['-g',groups]

# convert jser to Series object
in_series = Series.openJser(temp_jser)

# extract trace names
trace_names = in_series.objects.getNames()

def get_object_type(input_string, object_types=["sp", "ax", "c", "endo", "r", "ser", "sa"]):
    if 'endo' in input_string:
        input_string = input_string.split('_')[0]
    or_object_type_pattern = "|".join(object_types)
    pattern = fr'^(d([0-9]{{2}}))(({or_object_type_pattern})(h|n|b|s)?([0-9]{{2,3}}))$'
    match = re.match(pattern, input_string)
    
    # if we have an object that is associated with a dendrite
    if match is not None:
        return match.groups()[3], match.groups()[1]
        
    # otherwise we might have a dendrite
    pattern = r'^(d([0-9]{2}))$'
    match = re.match(pattern, input_string)
    if match is not None:
        return "d", match.groups()[1]
    
    # otherwise we don't really know what this object is.
    return None, None


group_names = {'sp':'spines','ax':'axons','d':'dendrites','c':'synapses','r':'ribosomes','ser':'ser','sa':'spine_apparatus','endo':'endosomes'}
neuron_dict = defaultdict(list)
hash_to_trace = {}
for trace in trace_names:
    obj,neuron = get_object_type(trace,list(group_names.keys()))
    if obj is not None:
        id = hashName(trace)
        hash_to_trace[id] = trace
        if obj == 'sp' or obj == 'd':
            neuron_dict[neuron].append(id)
        in_series.object_groups.add(group_names[obj],trace)

in_series.save()
in_series.saveJser()

in_series.close()

with open("lookup.json","w") as f:
    json.dump(neuron_dict, f, indent=4)

subprocess.run(run_command)

if 'neurons' in groups:
    # add neuron labels
    z = zarr.open(args.output)
    z.attrs['neuron_dict'] = neuron_dict

    d_arr = z['labels_dendrites'][:]
    s_arr = z["labels_spines"][:]

    neuron_labels = np.zeros_like(s_arr)

    for i in neuron_dict.keys():
        neuron_label = np.zeros_like(s_arr)
        ids = neuron_dict[i]
        for id in ids:
            mask = np.logical_or(d_arr==id, s_arr==id)
            neuron_labels[mask] = int(i)

    z['labels_neurons'] = neuron_labels
    z['labels_neurons'].attrs['resolution'] = z['labels_dendrites'].attrs['resolution']
    z['labels_neurons'].attrs['offset'] = z['labels_dendrites'].attrs['offset']