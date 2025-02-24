"""
Получение информации о варианте из ClinVar, используя esummary.

Скрипт получает информацию о генетическом варианте из базы данных ClinVar
и выводит эту информацию в консоль.

Скрипт использует API NCBI (esummary) и XML формат для обмена данными.

Использование:
    python clinvar_variant_info.py <variation_id>

Аргументы:
    <variation_id>  ClinVar Variation ID варианта, информацию о котором нужно получить (например, 9).

Пример:
    python clinvar_variant_info.py 9
    python clinvar_variant_info.py 1234

Автор: Smolin Igor
"""

import requests
import xml.etree.ElementTree as ET
import sys

def get_clinvar_info(variation_id):
    try:
        esummary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        esummary_params = {
            "db": "clinvar",
            "id": variation_id,
            "retmode": "xml"
        }
        esummary_response = requests.get(esummary_url, params=esummary_params)
        esummary_response.raise_for_status()  # Проверка на ошибки HTTP

        xml_content = esummary_response.content.decode('utf-8')
        root = ET.fromstring(xml_content)

        variant_info = {}

        doc_summary = root.find(".//DocumentSummary")
        if doc_summary is not None:
            variant_info["obj_type"] = doc_summary.findtext("obj_type")
            variant_info["accession"] = doc_summary.findtext("accession")
            variant_info["accession_version"] = doc_summary.findtext("accession_version")
            variant_info["title"] = doc_summary.findtext("title")

            germline_classification = doc_summary.find(".//germline_classification/description")
            if germline_classification is not None:
                variant_info["clinical_significance"] = germline_classification.text
                variant_info["clinical_significance_review_status"] = doc_summary.findtext(".//germline_classification/review_status") #Добавим статус проверки клинической значимости
            else:
                variant_info["clinical_significance"] = "Не найдено"
                variant_info["clinical_significance_review_status"] = "Не найдено"

            gene = doc_summary.find(".//genes/gene")
            if gene is not None:
                variant_info["gene_symbol"] = gene.findtext("symbol")
                variant_info["gene_id"] = gene.findtext("GeneID")
            else:
                variant_info["gene_symbol"] = "Не найдено"
                variant_info["gene_id"] = "Не найдено"

            scv_list = []
            for scv_element in doc_summary.findall(".//supporting_submissions/scv/string"):
                scv_list.append(scv_element.text)
            variant_info["supporting_submissions"] = scv_list

            rcv_list = []
            for rcv_element in doc_summary.findall(".//supporting_submissions/rcv/string"):
                rcv_list.append(rcv_element.text)
            variant_info["related_rcv_accessions"] = rcv_list #Более понятное название

            allele_frequencies = []
            for allele_freq_element in doc_summary.findall(".//allele_freq_set/allele_freq"):
                source = allele_freq_element.findtext("source")
                value = allele_freq_element.findtext("value")
                minor_allele = allele_freq_element.findtext("minor_allele")
                allele_frequencies.append({
                    "source": source,
                    "value": value,
                    "minor_allele": minor_allele
                })
            variant_info["allele_frequencies"] = allele_frequencies

            variation = doc_summary.find(".//variation_set/variation")
            if variation is not None:
                 variant_info["canonical_spdi"] = variation.findtext("canonical_spdi")
            else:
                variant_info["canonical_spdi"] = "Не найдено"

            if variation is not None:
                 variant_info["variant_type"] = variation.findtext("variant_type")
            else:
                variant_info["variant_type"] = "Не найдено"

            if variation is not None:
                 variant_info["cdna_change"] = variation.findtext("cdna_change")
            else:
                variant_info["cdna_change"] = "Не найдено"


            trait_set = []
            for trait_element in doc_summary.findall(".//germline_classification/trait_set/trait"):
                trait_name = trait_element.findtext("trait_name")
                trait_set.append(trait_name)
            variant_info["trait_set"] = trait_set

            return variant_info
        else:
            return "DocumentSummary не найден"

    except requests.exceptions.RequestException as e:
        return f"Ошибка запроса: {e}"
    except ET.ParseError as e:
        return f"Ошибка парсинга XML: {e}"
    except Exception as e:
        return f"Общая ошибка: {e}"

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Использование: python clinvar_variant_info.py <variation_id>")
        sys.exit(1)

    variation_id = sys.argv[1]  # Получаем variation_id из аргументов командной строки
    variant_data = get_clinvar_info(variation_id)

    if isinstance(variant_data, dict):
        print(f"\n----- Информация о варианте {variation_id} -----")
        print(f"Тип объекта: {variant_data.get('obj_type', 'Не найдено')}")
        print(f"Accession: {variant_data.get('accession', 'Не найдено')}")
        print(f"Accession версия: {variant_data.get('accession_version', 'Не найдено')}")
        print(f"Название: {variant_data.get('title', 'Не найдено')}")

        print(f"Клиническая значимость: {variant_data.get('clinical_significance', 'Не найдено')}")
        print(f"Статус проверки клинической значимости: {variant_data.get('clinical_significance_review_status', 'Не найдено')}")

        print(f"Символ гена: {variant_data.get('gene_symbol', 'Не найдено')}")
        print(f"Gene ID: {variant_data.get('gene_id', 'Не найдено')}")

        print(f"Список SCV: {variant_data.get('supporting_submissions', [])}")
        print(f"Список RCV: {variant_data.get('related_rcv_accessions', [])}")

        print(f"Canonical SPDI: {variant_data.get('canonical_spdi', 'Не найдено')}")
        print(f"Тип варианта: {variant_data.get('variant_type', 'Не найдено')}")
        print(f"cdna изменение: {variant_data.get('cdna_change', 'Не найдено')}")

        print("\nАллельные частоты:")
        for freq in variant_data.get('allele_frequencies', []):
            print(f"  Источник: {freq.get('source', 'Не найдено')}, Значение: {freq.get('value', 'Не найдено')}, Minor Allele: {freq.get('minor_allele', 'Не найдено')}")

        print("\nСписок связанных признаков/заболеваний:")
        for trait in variant_data.get('trait_set', []):
            print(f"  - {trait}")

        print("----- Конец информации о варианте -----")
    else:
        print(f"Ошибка при получении информации о варианте {variation_id}: {variant_data}")