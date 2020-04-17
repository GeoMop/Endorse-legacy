# projekt Endorse

Predikce vlastností EDZ (excavation damage zone) s vlivem na bezpečnost a spolehlivost hlubinného úložiště radioaktivního odpadu.

Obsah repozitáře:

- `mlmc-modelling-2019` - experimentální aplikace metody MLMC na puklinové modely, 2D proudění na čtverci
- `repository-model` - The repository scale model. Captures transport of a tracer along a single storage horizontal well
  with random large fractures, impact of the neighboring wells is considered.



## Cíl projektu

Vytvořit SW nástroj a metodiku, pro predikci veličin charakterizujících bezpečnost dílčí části úložiště
(tzv. *indikátorů bezpečnosti*) na základě geofyzikálních měření. To zahrnuje:

1. Sestavení modelu transportu kontaminace skrze EDZ od (náhodných) úložných kontejnerů do hypotetické poruchy. 
Zahrnutí předpokládané geometrie úložiště s velikostí do 100m.
2. Definice vhodných indikátorů bezpečnosti jakožto veličin odvozených od výzledků modelu transportu.
3. Tvorbu menších modelů pro identifikaci parametrů transportního modelu na základě předpokládaných průzkumů 
a geofyzikálních měření.
4. Aplikaci vhodných stochastických výpočetních metod pro predikci rozdělení indikátorů bezpečnosti a parametrů 
transportního modelu se zahrnutím chyb měření a dalších podstatných neurčitostí použitých modelů

## Rozcestník

- [Přehled řešení projektu](https://github.com/jbrezmorf/Endorse/projects/2) - přehled plánovaných, řešených a ukončených úkolů dle harmonogramu projektu

- [Přehled řešitelů](https://docs.google.com/document/d/1R8CBU9197brrruWGahVbE7_At2S2V51J6JV5bgs-kxQ/edit#heading=h.e1t1yg8nyvaz)

- [Zotero Endorse](https://www.zotero.org/groups/287302/flow123d/items/collectionKey/3BAS5Z2A) - sdílený prostor pro komunikaci referencí a fulltextů, použití v rámci aplikace [Zotero](https://www.zotero.org/download/)

- [Overleaf Endorse](https://www.overleaf.com/project) - tvorba sdílených textů, zpráv, ... 

## Software

- [Flow123d](https://github.com/flow123d/flow123d) 
 simulátor transportních a mechanických procesů v rozpukaném porézním prostředí

- [MLMC](https://github.com/GeoMop/MLMC)
  metoda multilevel Monte Carlo v Pythonu, generování náhodných polí a puklinových sítí, 
  maximal entropy method pro rekonstrukci hustoty pravděpodobnosti
  
- [PERMON](https://github.com/permon)
