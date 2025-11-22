# Template Structural Tags

Based on the [structural tags](./structural_tag.md), template structural tags provides a more convineient method for users to generate grammars to describe the grammar to constrain the LLMs' output.
In general, Template structural tags support placeholders in the structural tags, which are in the form of `{lists[].arg}`. These values will be automatically expanded with the user's input values. Users can use `xgrammar.Grammar.apply_structural_tag_template(template_json_str: Union[str, Dict[str, Any]], **kwargs: List[Dict[str, Any]])` to generate a `xgrammar.Grammar`. Here, `template_json_str` is the template structural tags, and kwargs are a series of values. For example, the template structural tag `stag` contains `{A[].name}` and `{A[].age}`, then the users can use `xgrammar.Grammar.apply_structural_tag_template(stag, A=A)` to generate the grammar, where `A=[{"name":..., "age":...},.{"name":...,"age":...}, ...]`. This function will replace the placeholders automatically.

## Template Placeholders in Formats

1. `const_string`

Each `const_string` format can contain multiple placeholders, but they must be from the same value mapping. for example, `A=[{"begin": "It is", "end":"."}, {"begin": "Is it", "end": "?"}]`, and the template is:

```json
{
    "type": "const_string",
    "value": "{A[].begin} a dog{A[].end}"
}
```

Is allowed. And the output is constrained to `It is a dog.` and `Is it a dog?`. However, if the provided values are `Begin=[{"word": "It is"}, {"word": "Is it"}], End=[{"word": "."}, {"word": "?"}]`, and the template is

```json
{
    "type": "const_string",
    "value": "{Begin[].word} a dog{End[].word}"
}
```

cannot be compiled. Because the meaning is ambiguous, and we call these templates **mingled**, we cannot compile them. This template format will be expanded into a `const_string` format or an `or` format, or return nothing.

2. `grammar`, `json_schema`, `qwen_xml_parameter`, `regex`

If the template placeholder is in these formats' value, **only if** the value is exactly the placeholder. For example,

```json
{
    "type": "json_schemas",
    "json_schema": "{schemas[].schema}"
}
```

can be automatically replaced with the given `schemas`. However, this format will not be replaced:

```json
{
    "type": "json_schemas",
    "json_schema": {
        "type": "{schemas[].schema}"
        }
}
```

The same rule holds for the four formats. This template format will be expanded into a `json_schemas` format, or an `or` format, or return nothing.

3. `tag`

For a tag, it is allowed to contain a placeholder in the `begin` and `end` fields. It is also okay if the `content` field is also a template format. However, as the same as `const_string`, `begin` and `end` can contain multiple placeholders, but they must be from the same value mapping. Otherwise, the template is mingled and cannot be compiled. For example, this is a valid tag template:

```json
{
    "type": "tag",
    "begin": "<function={tools[].name}",
    "content": {
        "type": "json_schema",
        "json_schema": "{tools[].args}"
        },
    "end": "</function>"
}
```

This format will be expanded into a `tag` format or an `or` format, or a series of `tag` formats in `triggered_tags`, `tags_with_separator`, or return nothing.

4. `triggered_tags`, `tags_with_separator`

Template placeholders are not allowed in `triggers` and `separators`.

5. `sequence`, `or`, `any_text`

Template placeholders cannot be directly contained in these formats.
## Valid Template Structural Tags
Not all the template structural tags are valid. For example, the mingled formats mentioned above are not valid template structural tags. Besides, there are some other situations where the template structural tag cannot be compiled. For example:

```json
{
    "type": "sequence",
    "elements": [
            {
            "type": "const_string",
            "value": "{A[].value1}"
            },
            {
            "type": "const_string",
            "value": "{B[].value1}"
            },
            {
            "type": "const_string",
            "value": "{A[].value2}"
            },
            {
            "type": "const_string",
            "value": "{B[].value2}"
            },
        ]
}
```

cannot be compiled because we cannot analyze the meaning of the sequence. However,

```json
{
    "type": "sequence",
    "elements": [
            {
                "type": "const_string",
                "value": "{A[].value1}"
            },
            {
                "type": "const_string",
                "value": "{B[].value}"
            },
            {
                "type": "const_string",
                "value": "{A[].value2}"
            },
        ]
}
```

can be compiled. It will be interpreted as a pair of `A.value1`, `A.value2`, and an arbitrary `B.value`. Basically, all the invalid template structural tags are in similar situations: we cannot divide each template placeholder's scope properly. In the former invalid one, both `A` and `B`'s scopes are the `sequence` format. In the latter valid one, `A`'s scope is the `sequence` format, while `B`'s scope is the `const_string` format.

## Builtin Template Structural Tags

There are several builtin template structural tags, designed for different LLMs' tool-calling formats. Users can get the builtin template structural tags with `xgrammar.structural_tag.get_builtin_template_structural_tag(format_type: str)`. Currently, these `format_type`s are supported: `Llama`, `Kimi`, `Deepseek`, `Qwen_Coder`, `Qwen`, `Harmony`.
For `Llama`, `Kimi`, `Deepseek`, `Qwen_Coder`, `Qwen`, users need to provide a value `tools=[{"name":..., "parameters":...}, ...]` for the template structural tags. For `Harmony` format, users must provide both `tools=[{"name":..., "parameters":...}]` and `builtin_tools=[{"name":..., "parameters":...}, ...]` for the template. These templates will force the LLMs to output the correct function-calling formats, with other natural language outputs.
