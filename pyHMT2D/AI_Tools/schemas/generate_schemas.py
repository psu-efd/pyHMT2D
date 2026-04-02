# -*- coding: utf-8 -*-
"""
Generate openai_tool_schemas.json from the AI_Tools functions.

Run this script once after any tool signature change:
    python -m pyHMT2D.AI_Tools.schemas.generate_schemas

The output file (openai_tool_schemas.json) is suitable for use with:
  - OpenAI function calling  (openai.ChatCompletion with `tools=` parameter)
  - Anthropic tool use       (anthropic.Anthropic with `tools=` parameter)
  - Any other LLM API that accepts JSON-Schema tool definitions
"""

from __future__ import annotations

import inspect
import json
import os
import typing
from typing import get_type_hints

import pyHMT2D.AI_Tools.tools as _tools


_TYPE_MAP = {
    "str":   {"type": "string"},
    "int":   {"type": "integer"},
    "float": {"type": "number"},
    "bool":  {"type": "boolean"},
    "list":  {"type": "array",  "items": {}},
    "dict":  {"type": "object"},
    "None":  {"type": "null"},
}


def _py_type_to_json_schema(annotation) -> dict:
    """Best-effort conversion of a Python type annotation to JSON Schema."""
    if annotation is inspect.Parameter.empty:
        return {}

    # Unwrap Optional[X] → X (since tools use Optional for optional params)
    origin = getattr(annotation, "__origin__", None)
    args   = getattr(annotation, "__args__", ())

    if origin is typing.Union and type(None) in args:
        # Optional[X]
        inner = [a for a in args if a is not type(None)]
        if len(inner) == 1:
            return _py_type_to_json_schema(inner[0])

    if origin is list or annotation is list:
        item_schema = {}
        if args:
            item_schema = _py_type_to_json_schema(args[0])
        return {"type": "array", "items": item_schema}

    if origin is dict or annotation is dict:
        return {"type": "object"}

    name = getattr(annotation, "__name__", str(annotation))
    return _TYPE_MAP.get(name, {"type": "string"})


def _build_tool_schema(fn) -> dict:
    """Build an OpenAI-compatible tool schema from a function."""
    sig    = inspect.signature(fn)
    hints  = get_type_hints(fn)
    doc    = inspect.getdoc(fn) or ""

    # Description: first paragraph of docstring
    description = doc.split("\n\n")[0].replace("\n", " ").strip()

    properties = {}
    required   = []

    for name, param in sig.parameters.items():
        if name in ("self", "cls"):
            continue

        annotation = hints.get(name, inspect.Parameter.empty)
        schema     = _py_type_to_json_schema(annotation)

        # Extract per-param description from docstring "Parameters" section
        param_desc = ""
        in_params  = False
        for line in doc.splitlines():
            stripped = line.strip()
            if stripped.lower().startswith("parameters"):
                in_params = True
                continue
            if in_params:
                if stripped.startswith("----------"):
                    continue
                if stripped.startswith(f"{name} :") or stripped.startswith(f"{name}:"):
                    # next non-empty line is the description
                    continue
                if stripped and not stripped.startswith("----"):
                    # Could be the description line or another param
                    if any(
                        stripped.startswith(p + " :") or stripped.startswith(p + ":")
                        for p in sig.parameters
                        if p != name
                    ):
                        break
                    if param_desc == "":
                        param_desc = stripped

        schema["description"] = param_desc or name

        if param.default is inspect.Parameter.empty:
            required.append(name)
        else:
            schema["default"] = (
                param.default if param.default is not None else None
            )

        properties[name] = schema

    return {
        "type": "function",
        "function": {
            "name": fn.__name__,
            "description": description,
            "parameters": {
                "type": "object",
                "properties": properties,
                "required": required,
            },
        },
    }


def generate() -> None:
    from pyHMT2D.AI_Tools.tools import __all__ as tool_names
    import pyHMT2D.AI_Tools.tools as tool_module

    schemas = []
    for name in tool_names:
        fn = getattr(tool_module, name, None)
        if fn and callable(fn):
            schemas.append(_build_tool_schema(fn))

    out_path = os.path.join(os.path.dirname(__file__), "openai_tool_schemas.json")
    with open(out_path, "w") as fh:
        json.dump(schemas, fh, indent=2)

    print(f"Generated {len(schemas)} tool schemas → {out_path}")


if __name__ == "__main__":
    generate()
