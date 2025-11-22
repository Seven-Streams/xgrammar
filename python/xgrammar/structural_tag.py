"""Defines all structural tag formats."""

import json
from typing import Any, Dict, List, Literal, Type, Union

try:
    # Python 3.9+
    from typing import Annotated
except ImportError:
    # Python 3.8
    from typing_extensions import Annotated

from pydantic import BaseModel, Field

# ---------- Basic Formats ----------


class ConstStringFormat(BaseModel):
    """A format that matches a constant string."""

    type: Literal["const_string"] = "const_string"
    """The type of the format."""
    value: str
    """The constant string."""


class JSONSchemaFormat(BaseModel):
    """A format that matches a JSON schema."""

    type: Literal["json_schema"] = "json_schema"
    """The type of the format."""
    json_schema: Union[bool, Dict[str, Any]]
    """The JSON schema."""


class QwenXMLParameterFormat(BaseModel):
    """A format that matches Qwen XML function calls.

    Examples
    --------
    .. code-block:: python

        structural_tag = QwenXMLParameterFormat(
            json_schema={
                "type": "qwen_xml_parameter",
                "json_schema": {
                    "type": "object",
                    "properties": {"name": {"type": "string"}, "age": {"type": "integer"}},
                    "required": ["name", "age"],
                },
            }
        )

    The above structural tag can accept the following outputs::

        <parameter=name>Bob</parameter><parameter=age>100</parameter>
        <parameter=name>"Bob&lt;"</parameter><parameter=age>100</parameter>

    """

    type: Literal["qwen_xml_parameter"] = "qwen_xml_parameter"
    """The type of the format."""

    json_schema: Union[bool, Dict[str, Any]]
    """The JSON schema for the parameters of the function calling."""


class AnyTextFormat(BaseModel):
    """A format that matches any text."""

    type: Literal["any_text"] = "any_text"
    """The type of the format."""


class GrammarFormat(BaseModel):
    """A format that matches an ebnf grammar."""

    type: Literal["grammar"] = "grammar"
    """The type of the format."""

    grammar: str
    """The ebnf grammar."""


class RegexFormat(BaseModel):
    """A format that matches a regex pattern."""

    type: Literal["regex"] = "regex"
    """The type of the format."""

    pattern: str
    """The regex pattern."""


# ---------- Combinatorial Formats ----------


class SequenceFormat(BaseModel):
    """A format that matches a sequence of formats."""

    type: Literal["sequence"] = "sequence"
    """The type of the format."""
    elements: List["Format"]
    """The elements of the sequence."""


class OrFormat(BaseModel):
    """A format that matches one of the formats."""

    type: Literal["or"] = "or"
    """The type of the format."""
    elements: List["Format"]
    """The elements of the or."""


class TagFormat(BaseModel):
    """A format that matches a tag: ``begin content end``."""

    type: Literal["tag"] = "tag"
    """The type of the format."""
    begin: str
    """The begin tag."""
    content: "Format"
    """The content of the tag. It can be any of the formats."""
    end: str
    """The end tag."""


class TriggeredTagsFormat(BaseModel):
    """A format that matches triggered tags. It can allow any output until a trigger is
    encountered, then dispatch to the corresponding tag; when the end tag is encountered, the
    grammar will allow any following output, until the next trigger is encountered.

    Each tag should be matched by exactly one trigger. "matching" means the trigger should be a
    prefix of the begin tag.

    Examples
    --------

    .. code-block:: python

        structural_tag = TriggeredTagsFormat(
            triggers=["<function="],
            tags=[
                TagFormat(
                    begin="<function=func1>",
                    content=JSONSchemaFormat(json_schema=...),
                    end="</function>",
                ),
                TagFormat(
                    begin="<function=func2>",
                    content=JSONSchemaFormat(json_schema=...),
                    end="</function>",
                ),
            ],
            at_least_one=False,
            stop_after_first=False,
        )

    The above structural tag can accept the following outputs::

        <function=func1>{"name": "John", "age": 30}</function>
        <function=func2>{"name": "Jane", "age": 25}</function>
        any_text<function=func1>{"name": "John", "age": 30}</function>any_text1<function=func2>{"name": "Jane", "age": 25}</function>any_text2

    """

    type: Literal["triggered_tags"] = "triggered_tags"
    """The type of the format."""
    triggers: List[str]
    """The triggers of the triggered tags."""
    tags: List[TagFormat]
    """The tags of the triggered tags."""
    at_least_one: bool = False
    """Whether at least one of the tags must be generated."""
    stop_after_first: bool = False
    """Whether to stop after the first tag is generated."""


class TagsWithSeparatorFormat(BaseModel):
    """A format that matches a tags with separator. It can match zero, one, or more tags, separated
    by the separator, with no other text allowed.

    Examples
    --------

    .. code-block:: python

        structural_tag = TagsWithSeparatorFormat(
            tags=[
                TagFormat(begin="<function=func1>", content=JSONSchemaFormat(json_schema=...), end="</function>"),
                TagFormat(begin="<function=func2>", content=JSONSchemaFormat(json_schema=...), end="</function>"),
            ],
            separator=",",
            at_least_one=False,
            stop_after_first=False,
        )

    The above structural tag can accept an empty string, or the following outputs::

        <function=func1>{"name": "John", "age": 30}</function>
        <function=func1>{"name": "John", "age": 30}</function>,<function=func2>{"name": "Jane", "age": 25}</function>
        <function=func1>{"name": "John", "age": 30}</function>,<function=func2>{"name": "Jane", "age": 25}</function>,<function=func1>{"name": "John", "age": 30}</function>
    """

    type: Literal["tags_with_separator"] = "tags_with_separator"
    """The type of the format."""
    tags: List[TagFormat]
    """The tags of the tags with separator."""
    separator: str
    """The separator of the tags with separator."""
    at_least_one: bool = False
    """Whether at least one of the tags must be matched."""
    stop_after_first: bool = False
    """Whether to stop after the first tag is matched."""


# ---------- Discriminated Union ----------


Format = Annotated[
    Union[
        AnyTextFormat,
        ConstStringFormat,
        JSONSchemaFormat,
        GrammarFormat,
        RegexFormat,
        QwenXMLParameterFormat,
        OrFormat,
        SequenceFormat,
        TagFormat,
        TriggeredTagsFormat,
        TagsWithSeparatorFormat,
    ],
    Field(discriminator="type"),
]
"""Union of all structural tag formats."""


# Solve forward references
if hasattr(BaseModel, "model_rebuild"):
    SequenceFormat.model_rebuild()
    TagFormat.model_rebuild()
    TriggeredTagsFormat.model_rebuild()
    TagsWithSeparatorFormat.model_rebuild()
elif hasattr(BaseModel, "update_forward_refs"):
    # This is for backward compatibility with pydantic v1
    SequenceFormat.update_forward_refs()
    TagFormat.update_forward_refs()
    TriggeredTagsFormat.update_forward_refs()
    TagsWithSeparatorFormat.update_forward_refs()
else:
    raise RuntimeError("Unsupported pydantic version")


# ---------- Top Level ----------


class StructuralTagItem(BaseModel):
    """Deprecated. Definition of a structural tag item.

    See :meth:`xgrammar.Grammar.from_structural_tag` for more details.
    """

    begin: str
    """The begin tag."""
    schema_: Union[str, Type[BaseModel], Dict[str, Any]] = Field(alias="schema")
    """The schema."""
    end: str
    """The end tag."""


class StructuralTag(BaseModel):
    """
    Describes a complete structural tag structure. It corresponds to
    ``"response_format": {"type": "structural_tag", "format": {...}}`` in API.
    """

    type: Literal["structural_tag"] = "structural_tag"
    """The type must be "structural_tag"."""
    format: Format
    """The format of the structural tag. Could be any of the structural tag formats."""

    @staticmethod
    def from_legacy_structural_tag(
        tags: List[StructuralTagItem], triggers: List[str]
    ) -> "StructuralTag":
        """Convert a legacy structural tag item to a structural tag."""
        return StructuralTag(
            type="structural_tag",
            format=TriggeredTagsFormat(
                type="triggered_tags",
                triggers=triggers,
                tags=[
                    TagFormat(
                        begin=tag.begin,
                        content=JSONSchemaFormat(
                            json_schema=(
                                json.loads(tag.schema_)
                                if isinstance(tag.schema_, str)
                                else (
                                    tag.schema_.model_json_schema()
                                    if isinstance(tag.schema_, type)
                                    and issubclass(tag.schema_, BaseModel)
                                    else tag.schema_
                                )
                            )
                        ),
                        end=tag.end,
                    )
                    for tag in tags
                ],
            ),
        )

    @staticmethod
    def from_json(json_str: Union[str, Dict[str, Any]]) -> "StructuralTag":
        """Convert a JSON string to a structural tag."""
        if isinstance(json_str, str):
            return StructuralTag.model_validate_json(json_str)
        elif isinstance(json_str, dict):
            return StructuralTag.model_validate(json_str)
        else:
            raise ValueError("Invalid JSON string or dictionary")


__all__ = [
    "ConstStringFormat",
    "JSONSchemaFormat",
    "QwenXMLParameterFormat",
    "AnyTextFormat",
    "GrammarFormat",
    "RegexFormat",
    "SequenceFormat",
    "OrFormat",
    "TagFormat",
    "TriggeredTagsFormat",
    "TagsWithSeparatorFormat",
    "Format",
    "StructuralTagItem",
    "StructuralTag",
]

# ---------- Template Structural Tag ----------

_structural_tag_registry = {}


def _register_template_structural_tag_format(name: str):
    """Register a structural tag format."""

    def decorator(func):
        _structural_tag_registry[name] = func
        return func

    return decorator


def get_builtin_template_structural_tag(format_type: str) -> Dict:
    """Get builtin template structural tag format by format type.
    In all the template structural tag formats, users should provide
    a list of tools, each tool should have a "name" and "parameters" field.
    to use the template structural tag format. Besides, for the OpenAI Harmony Response Format,
    users should also provide a list of builtin tools, each builtin tool should have a "name"
    and "parameters" field.

    Examples
    --------

    .. code-block:: python

        from xgrammar import structural_tag, Grammar
        tools = [
            {"name": "tool1", "parameters": {"param1": {"type": "string"}}},
            {"name": "tool2", "parameters": {"param2": {"type": "integer"}}},
        ]
        builtin_tools = [
            {"name": "builtin_tool1", "parameters": {"param1": {"type": "string"}}},
            {"name": "builtin_tool2", "parameters": {"param2": {"type": "integer"}}},
        ]
        template_structural_tag = structural_tag.get_builtin_template_structural_tag("Harmony")
        grammar = Grammar.apply_template_structural_tag(template_structural_tag, tools=tools, builtin_tools=builtin_tools)

    The above grammar can be used to construct a grammar that matches the function calling
    format of the specified model.



    Parameters
    ----------
    format_type : str
        The format type, must be one of the registered format types:
        "Llama", "Kimi", "Deepseek", "Qwen_Coder", "Qwen", "Harmony".

    Returns
    -------
    Dict
        A template structural tag format dictionary.

    """
    func = _structural_tag_registry.get(format_type)
    if func is None:
        raise ValueError(f"Unknown format type: {format_type}")
    return func()


@_register_template_structural_tag_format("Llama")
def get_llama_style_template_structural_tag() -> Dict:
    """Get Llama style structural tag format.

    Returns
    -------
    Dict
        A template structural tag format dictionary.
        This format is used by Llama 3 and other models that follow the same style.
    """
    return {
        "type": "structural_tag",
        "format": {
            "type": "triggered_tags",
            "triggers": ['{"name": '],
            "tags": [
                {
                    "begin": '{"name": "{tools[].name}", "parameters": ',
                    "content": {"type": "json_schema", "json_schema": "{tools[].parameters}"},
                    "end": "}",
                }
            ],
        },
    }


@_register_template_structural_tag_format("Kimi")
def get_kimi_style_template_structural_tag() -> Dict:
    """Get Kimi style structural tag format.

    Returns
    -------
    Dict
        A template structural tag format dictionary.
        This format is used by Kimi-v2 and other models that follow the same style.
    """
    return {
        "type": "structural_tag",
        "format": {
            "type": "triggered_tags",
            "triggers": ["<|tool_call_begin|>"],
            "tags": [
                {
                    "begin": "<|tool_call_begin|>{tools[].name}<|tool_call_argument_begin|>",
                    "content": {"type": "json_schema", "json_schema": "{tools[].parameters}"},
                    "end": "<|tool_call_end|>",
                }
            ],
        },
    }


@_register_template_structural_tag_format("Deepseek")
def get_deepseek_style_template_structural_tag() -> Dict:
    """Get Deepseek style structural tag format.

    Returns
    -------
    Dict
        A template structural tag format dictionary.
        This format is used by Deepseek-v3.1 and other models that follow the same style.
    """
    return {
        "type": "structural_tag",
        "format": {
            "type": "triggered_tags",
            "triggers": ["<｜tool▁calls▁begin｜><｜tool▁call▁begin｜>"],
            "tags": [
                {
                    "begin": "<｜tool▁calls▁begin｜><｜tool▁call▁begin｜>{tools[].name}<｜tool▁sep｜>",
                    "content": {"type": "json_schema", "json_schema": "{tools[].parameters}"},
                    "end": "<｜tool▁call▁end｜>",
                }
            ],
        },
    }


@_register_template_structural_tag_format("Qwen_Coder")
def get_qwen_coder_style_template_structural_tag() -> Dict:
    """Get Qwen Coder style structural tag format.

    Returns
    -------
    Dict
        A template structural tag format dictionary.
        This format is used by Qwen3 Coder and other models that follow the same style.
    """
    return {
        "type": "structural_tag",
        "format": {
            "type": "triggered_tags",
            "triggers": ["<function="],
            "tags": [
                {
                    "begin": "<function={tools[].name}>",
                    "content": {
                        "type": "qwen_xml_parameter",
                        "json_schema": "{tools[].parameters}",
                    },
                    "end": "</function>",
                }
            ],
        },
    }


@_register_template_structural_tag_format("Qwen")
def get_qwen_style_template_structural_tag() -> Dict:
    """Get Qwen style structural tag format.

    Returns
    -------
    Dict
        A template structural tag format dictionary.
        This format is used by Qwen3 and other models that follow the same style.
    """
    return {
        "type": "structural_tag",
        "format": {
            "type": "triggered_tags",
            "triggers": ["<tool_call>"],
            "tags": [
                {
                    "begin": '<tool_call>{"name": "{tools[].name}", "arguments": ',
                    "content": {"type": "json_schema", "json_schema": "{tools[].parameters}"},
                    "end": "}</tool_call>",
                }
            ],
        },
    }


@_register_template_structural_tag_format("Harmony")
def get_harmony_style_template_structural_tag() -> Dict:
    """Get harmony style structural tag format.

    Returns
    -------
    Dict
        A template structural tag format dictionary.
        This format is in OpenAI Harmony Response Format, which is used by GPT-oss
        and other models that follow the same style.
    """
    return {
        "type": "structural_tag",
        "format": {
            "type": "triggered_tags",
            "triggers": ["<|start|>"],
            "tags": [
                {
                    "type": "tag",
                    "begin": "<|start|>assistant<|channel|>analysis<|message|>",
                    "content": {"type": "any_text"},
                    "end": "<|end|>",
                },
                {
                    "type": "tag",
                    "begin": "<|start|>assistant<|channel|>final<|message|>",
                    "content": {"type": "any_text"},
                    "end": "<|return|>",
                },
                {
                    "type": "tag",
                    "begin": "<|start|>assistant<|channel|>final<|message|>",
                    "content": {"type": "any_text"},
                    "end": "<|call|>",
                },
                {
                    "type": "tag",
                    "begin": "<|start|>assistant<|channel|>commentary to={tools[].name}<|constrain|>json<|message|>",
                    "content": {"type": "json_schema", "json_schema": "{tools[].parameters}"},
                    "end": "<|end|>",
                },
                {
                    "type": "tag",
                    "begin": "<|start|>assistant<|channel|>analysis to={builtin_tools[].name}<|message|>",
                    "content": {
                        "type": "json_schema",
                        "json_schema": "{builtin_tools[].parameters}",
                    },
                    "end": "<|end|>",
                },
            ],
        },
    }
