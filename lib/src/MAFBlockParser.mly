%{
%}
%token <string> A
%token <MAFDefs.sequence> S
%token BLANK EOF
%start parse parse_minimally
%type <MAFDefs.block> parse
%type <MAFDefs.block> parse_minimally
%%
parse:
    | block BLANK           { $1 }
    | block EOF             { $1 }
    | BLANK parse           { $2 }
;
block:
    | A sequences           { { MAFDefs.attributes = KeyValueParser.parse KeyValueLexer.token
                                                        (Lexing.from_string $1);
                                sequences = Array.of_list (List.rev $2);
                                unparsed = "" } }
;
sequences:
    | S                     { [$1] }
    | sequences S           { $2 :: $1 }
;
parse_minimally:
    | A S                   { { MAFDefs.attributes = [];
                                sequences = [| $2 |];
                                unparsed = "" } }
%%
