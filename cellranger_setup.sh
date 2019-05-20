# download cellranger software
curl -o cellranger-3.0.2.tar.gz "http://cf.10xgenomics.com/releases/cell-exp/cellranger-3.0.2.tar.gz?Expires=1558433017&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cDovL2NmLjEweGdlbm9taWNzLmNvbS9yZWxlYXNlcy9jZWxsLWV4cC9jZWxscmFuZ2VyLTMuMC4yLnRhci5neiIsIkNvbmRpdGlvbiI6eyJEYXRlTGVzc1RoYW4iOnsiQVdTOkVwb2NoVGltZSI6MTU1ODQzMzAxN319fV19&Signature=hBz4oscfE0mMJr49Aa2Aw00KBL-k094we4EznC2r~VNf~8fAUwz7ZFxeKGdUSqKtyVn-kMZxz6uQk4v1dUO2fmooeoS6ApN8DAcRESu6s69OSJVQyhv8I7zwSB5yZ7yGXLm~B3FYXRPp3tgrPIPmptGJ1rREZ2uwjdGZ5XkY2518iLIyB4AQlpg6U9ZRZvIvlp4lGsZ1gBE9-dGihmoDFTbfENEkmhvcGL8TK8d92KRzSA3nOCcBxchxAIKlXfWFDCQaehP658HJDi6xHcXo6G-2~bQqmmaj~1k3Kl~~yKP-wgR82wHGeiKxHkvtEPYkxp-57e0krrNQ6lA4cqLPzA__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
tar -xzvf cellranger-3.0.2.tar.gz

# download references
curl -O http://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-GRCh38-3.0.0.tar.gz

curl -O http://cf.10xgenomics.com/supp/cell-vdj/refdata-cellranger-vdj-GRCh38-alts-ensembl-2.0.0.tar.gz
tar -xzvf refdata-cellranger-vdj-GRCh38-alts-ensembl-2.0.0.tar.gz
