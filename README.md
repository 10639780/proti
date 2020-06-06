# proti
Proteine folding tool

to-do: inlezen
Bron over de rol van tijd in eiwitsynthese: met name het kopje 'Quality-Control Mechanisms Operate at Many Stages of Translation'

https://www.ncbi.nlm.nih.gov/books/NBK26829/#:~:text=The%20synthesis%20of%20most%20protein,each%20mRNA%20molecule%20being%20translated.

to-do: constraints en state-space begrijpen (college's van daan toepassen op case).

De state space van het kleinste molecuul is 3^8= 6561. Omdat het op de eerste plek waar gevouwen kan worden arbitrair is of er linksaf danwel rechtsaf gevouwen wordt is de state space direct afgenomen met 1/3. Hiermee nemen we aan dat spiegelbeelden van dezelfde moleculen als gelijkwaardig worden beschouwd (NB in de werkelijkheid kunnen spiegelbeeld isomeren zich wel degelijk anders gedragen). In onze benadering van de case zullen we ons beperken tot de regels die de stabiliteit van het eiwit bepalen. H-bindingen en later C-C (-5) en C-H (-1) bindingen.

Als we dit inzicht als constraint aannemen kunnen we de berekenbaarheid van het probleem verbeteren zonder een mogelijk optimale uitkomst weg te gooien. Een ander voorbeeld is de hard constraint dat het molecuul zichzelf niet mag kruisen. Kunnen we benaderen hoeveel mogelijkheden hiermee worden weggegooid? 

