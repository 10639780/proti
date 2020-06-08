# proti
Proteine folding tool

Bron over de rol van tijd in eiwitsynthese: met name het kopje 'Quality-Control Mechanisms Operate at Many Stages of Translation'

https://www.ncbi.nlm.nih.gov/books/NBK26829/#:~:text=The%20synthesis%20of%20most%20protein,each%20mRNA%20molecule%20being%20translated.

to-do: constraints en state-space begrijpen (college's van daan en bas toepassen op case).

Als spiegelbeelden met een constraint wegfilteren kunnen we de berekenbaarheid van het probleem verbeteren zonder een mogelijk optimale uitkomst weg te gooien. 

Wat betreft de harde constraint dat het molecuul zichzelf niet mag kruisen, kunnen we benaderen hoeveel mogelijkheden hiermee worden weggegooid? 

monteko.py vindt bij lengte 8 met 500 iteraties, in 40 minuten: 194 unieke eiwitten met allemaal een optimale score van -6.
