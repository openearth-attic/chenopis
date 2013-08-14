import os
import json
import collections
import itertools

import mako
from mako.template import Template

# Default location of variable file
SRCDIR = "."
VARFILE = os.path.join(SRCDIR, "variables.json")

# Read json file
def read_variables(filename=VARFILE):
    """Read the variables json file, append defaults from the same file"""
    with open(filename) as f:
        data = json.load(f)
    defaults = data['default']
    vars = []
    for variable in data['variables']:
        var = {}
        var.update(defaults)
        var.update(variable)
        if not var['use']:
            continue
        vars.append(var)
    return vars
vars = read_variables()

def read_modules(filename=VARFILE):
    """Read the variables json file, append defaults from the same file"""
    with open(filename) as f:
        data = json.load(f)
        modules = data['modules']
        return set(modules)
extra_modules = set(itertools.chain(var["modules"] for var in vars if var['modules']))
modules = read_modules() | extra_modules

# Templates, just keep them here for now
rank_template=\
"""
select case(var_name)
% for var in vars:
  case("${var["name"]}")
    rank = ${len(var["shape"])}
% endfor
  case default
    rank = 0
end select
"""

shape_template=\
"""
select case(var_name)
% for var in vars:
  case("${var["name"]}")
% if (var["shape"]):
    shape(1:${len(var["shape"])}) = (/${", ".join(var["shape"])}/)
% else:
    shape(:) = 0
% endif
% endfor
  case default
    shape(:) = 0
end select
"""

type_template=\
"""
select case(var_name)
% for var in vars:
  case("${var["name"]}")
    type_name = "${var["type"]}"
% endfor
  case default
    type_name = ""
end select
"""


unit_template=\
               """
select case(var_name)
% for var in vars:
  case("${var["name"]}")
    unit = "${var["unit"]}"
% endfor
  case default
    unit = ""
end select
"""


module_template=\
"""
% for module in modules:
use ${module}
% endfor
"""


nvar_template=\
"""
n = ${len(vars)}
"""

var_name_template=\
"""
select case(i)
% for i, var in enumerate(vars):
case(${i})
  name = "${var["name"]}"
% endfor
case default
  name = ""
end select
"""

# Generate a list of pointers
get_nd_template=\
"""
select case(var_name)
% for var in (var for var in vars):
<% ndvar = "x_{}d_{}_ptr".format(len(var["shape"]), var["type"]) %>
<% internalvar = var.get("internal", var["name"]) %>
% if var["pointer"]:
case("${var["name"]}")
  x = c_loc(${internalvar})
% else:
case("${var["name"]}")
! Allocate for nd >= 1
% if len(var["shape"]) > 0:
  if (allocated(${ndvar})) deallocate(${ndvar})
<% shapeargument = ",".join("size({},{})".format(internalvar, i + 1) for i in range(len(var["shape"]))) %>
  allocate(${ndvar}(${shapeargument}))
% endif
  ${ndvar} =  ${internalvar}
  x = c_loc(${ndvar})
% endif
% endfor
end select
"""


# Generate a list of pointers
set_nd_template=\
"""
select case(var_name)
% for var in (var for var in vars):
<% ndvar = "x_{}d_{}_ptr".format(len(var["shape"]), var["type"]) %>
<% internalvar = var.get("internal", var["name"]) %>
<% index = "(" + ",".join(":"*len(var["shape"])) +")" if (":" not in internalvar and len(var["shape"]) > 0)  else "" %>
case("${var['name']}")
% if var["pointer"]:
  call c_f_pointer(xptr, ${ndvar}, shape(${internalvar}))
  ${internalvar}${index} = ${ndvar}
% else:
  call c_f_pointer(xptr, ${ndvar}, shape(${internalvar}))
  ${internalvar}${index} = ${ndvar}
% endif
% endfor
end select
"""

# Which templates do we have...
templates = {
    "rank": rank_template,
    "type": type_template,
    "shape": shape_template,
    "module": module_template,
    "get_nd": get_nd_template,
    "set_nd": set_nd_template,
    "unit": unit_template,
    "nvar": nvar_template,
    "var_name": var_name_template
}


for include in templates.keys():
    outfilename = os.path.join(SRCDIR, "bmi_{}.inc".format(include))
    with open(outfilename, "w") as outfile:
        rendered = Template(templates[include]).render(vars=vars, modules=modules)
        outfile.write(rendered)
        print("!"*50)
        print(include)
        print("!"*50)
        print(rendered)
