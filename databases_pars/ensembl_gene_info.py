"""
Получение информации о гене из Ensembl API по символу гена и виду.

Скрипт получает информацию о гене из базы данных Ensembl
и выводит эту информацию в консоль.

Использование:
    python ensembl_gene_info.py <вид> <символ_гена>

Аргументы:
    <вид>          Вид, для которого нужно получить информацию (например, homo_sapiens, mus_musculus).
    <символ_гена>  Символ гена, информацию о котором нужно получить (например, BRCA2, Trp53).

Пример:
    python ensembl_gene_info.py homo_sapiens TP53
    python ensembl_gene_info.py mus_musculus Trp53

Автор: Smolin Igor
"""

import requests
import json
import sys

def get_gene_info_ensembl(species, gene_symbol):
    server = "https://rest.ensembl.org"
    ext = f"/lookup/symbol/{species}/{gene_symbol}?expand=1"
    try:
        response = requests.get(server + ext, headers={"Content-Type": "application/json"})
        response.raise_for_status()  # Raise HTTPError for bad responses

        if response.status_code == 200:
            gene_info = response.json()
            return gene_info
        else:
            print(f"Ошибка: {response.status_code} - {response.text}")
            return None

    except requests.exceptions.RequestException as e:
        print(f"Ошибка при запросе к Ensembl API: {e}")
        return None
    except json.JSONDecodeError as e:
        print(f"Ошибка при разборе JSON: {e}")
        return None

def print_gene_info(gene_info):
    if gene_info:
        print(f"\n----- Информация о гене {gene_info.get('display_name', 'Не найдено')} -----")
        print(f"Ensembl ID: {gene_info.get('id', 'Не найдено')}")
        print(f"Описание: {gene_info.get('description', 'Не найдено')}")
        print(f"Хромосома: {gene_info.get('seq_region_name', 'Не найдено')}")
        print(f"Старт: {gene_info.get('start', 'Не найдено')}")
        print(f"Стоп: {gene_info.get('end', 'Не найдено')}")
        print(f"Направление: {'прямая' if gene_info.get('strand', 0) == 1 else 'обратная'}")
        print(f"Тип: {gene_info.get('biotype', 'Не найдено')}")
        print(f"Источник: {gene_info.get('source', 'Не найдено')}")

        print("\nТранскрипты:")
        if 'Transcript' in gene_info:
            for transcript in gene_info['Transcript']:
                print(f"  Транскрипт ID: {transcript.get('id', 'Не найдено')}")
                print("    Экзоны:")
                if 'Exon' in transcript:
                    for exon in transcript['Exon']:
                        print(f"      Exon ID: {exon.get('id', 'Не найдено')}")
                        print(f"        Старт: {exon.get('start', 'Не найдено')}")
                        print(f"        Стоп: {exon.get('end', 'Не найдено')}")
                        print(f"        Хромосома: {exon.get('seq_region_name', 'Не найдено')}")
                        print(f"        Направление: {'прямая' if exon.get('strand', 0) == 1 else 'обратная'}")
                else:
                    print("      Экзоны не найдены.")
        else:
            print("    Транскрипты не найдены.")

        print("----- Конец информации о гене -----")
    else:
        print("Информация о гене не найдена.")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Использование: python ensembl_gene_info.py <вид> <символ_гена>")
        sys.exit(1)

    species = sys.argv[1]
    gene_symbol = sys.argv[2]
    gene_info = get_gene_info_ensembl(species, gene_symbol)

    if gene_info:
       print_gene_info(gene_info)
    else:
        print(f"Не удалось получить информацию о гене {gene_symbol} для вида {species} из Ensembl.")