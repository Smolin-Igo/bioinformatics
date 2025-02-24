import matplotlib
matplotlib.use('QtAgg')
import matplotlib.pyplot as plt
import numpy as np
import requests
from scipy.stats import chisquare

def get_fasta_from_ncbi(accession_number, email="your_email@example.com"):
    """Получает последовательность FASTA из NCBI по Accession Number."""
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "nucleotide",
        "id": accession_number,
        "rettype": "fasta",
        "retmode": "text",
        "email": email  # Обязательно укажите свой email
    }
    try:
        response = requests.get(base_url, params=params)
        response.raise_for_status()  # Raises HTTPError for bad requests (4XX, 5XX)
        fasta_string = response.text
        return fasta_string
    except requests.exceptions.RequestException as e:
        print(f"Ошибка при получении последовательности для {accession_number}: {e}")
        return None

def read_fasta_string(fasta_string):
    """Читает FASTA формат из строки и возвращает словарь."""
    sequences = {}
    try:
        header = None
        sequence = ''
        for line in fasta_string.splitlines():  # Разбиваем строку на строки
            line = line.strip()
            if line.startswith('>'):
                if header:
                    sequences[header] = sequence
                header = line[1:]
                sequence = ''
            else:
                sequence += line
        if header:
            sequences[header] = sequence
        return sequences
    except Exception as e:  # Указываем конкретный тип исключения
        print(f"Ошибка при парсинге FASTA строки: {e}")
        return {}

def calculate_gc_content(sequence):
    """Рассчитывает GC-состав последовательности ДНК."""
    sequence = sequence.upper()
    gc_count = sequence.count('G') + sequence.count('C')
    length = len(sequence)
    if length == 0:
        return 0.0
    return (gc_count / length) * 100

def calculate_nucleotide_counts(sequence):
    """Рассчитывает количество каждого нуклеотида в последовательности."""
    sequence = sequence.upper()
    a_count = sequence.count('A')
    c_count = sequence.count('C')
    g_count = sequence.count('G')
    t_count = sequence.count('T')
    return a_count, c_count, g_count, t_count

def perform_chisquare_test(sequence, expected_gc=50.0):
    """Проверяет GC-состав на отклонение от ожидаемого значения и возвращает результаты, включая observed_gc."""
    a_count, c_count, g_count, t_count = calculate_nucleotide_counts(sequence)
    observed = [c_count, g_count, a_count, t_count]
    total = sum(observed)

    expected_gc_count = total * (expected_gc / 100) / 2
    expected_at_count = total * ((100 - expected_gc) / 100) / 2
    expected = [expected_gc_count, expected_gc_count, expected_at_count, expected_at_count]

    observed_filtered = [o for o, e in zip(observed, expected) if e != 0]
    expected_filtered = [e for o, e in zip(observed, expected) if e != 0]

    if len(observed_filtered) < 2:
        return None, None, None, None, None  # Недостаточно данных

    statistic, p_value = chisquare(f_obs=observed_filtered, f_exp=expected_filtered)

    observed_gc = (c_count + g_count) / total * 100 if total > 0 else 0.0
    return statistic, p_value, observed_gc, expected_gc, total  # Возвращаем и expected_gc, и длину последовательности

def visualize_gc_content(header, accession_number, observed_gc, expected_gc, statistic, p_value, length):
    """Визуализирует GC-состав с пометкой о статистической значимости и длиной последовательности."""
    labels = ['Реальный', 'Ожидаемый']
    gc_values = [observed_gc, expected_gc]
    x = np.arange(len(labels))
    width = 0.35

    fig, ax = plt.subplots()
    rects = ax.bar(x, gc_values, width, label='GC-состав')

    # Добавляем пометку о статистической значимости
    if statistic is not None and p_value is not None:
        significance_label = "Значимо" if p_value < 0.05 else "Не значимо"
        ax.text(0.5, 0.85, f"Различия: {significance_label}", ha='center', va='top', transform=ax.transAxes, color='red') #Немного опустил для красоты

    # Добавляем пометку о длине последовательности
    ax.text(0.5, 0.95, f"Длина: {length}", ha='center', va='top', transform=ax.transAxes, color='blue')

    ax.set_ylabel('GC-состав (%)')
    ax.set_title(f'GC-состав: {header} ({accession_number})')
    ax.set_xticks(x, labels)
    ax.legend()
    ax.bar_label(rects, padding=3, fmt='%.2f')

    fig.tight_layout()
    plt.savefig(f"{accession_number}_gc_content.png")  # Сохраняем график
    plt.draw()  # Показываем график
    plt.pause(2)  # Пауза на 2 секунды
    plt.close(fig)  # Закрываем график


def main():
    """Основная функция для анализа GC-состава."""
    accession_numbers = [
        "NM_000518.5",  # Human TP53 mRNA
        "NC_000913.3",  # E. coli genome
        "NM_001363742.1", # Mouse Xist lncRNA
        "NC_007346.1",  # HIV-1 genome
        "NM_001304403.2", # Zebrafish actin beta 2, mRNA
        "NC_001422.1",  # Bacteriophage lambda, complete genome
        "NR_046018.1"   # Human microRNA 122-5p
    ]

    print("Анализ GC-состава последовательностей ДНК из NCBI\n")

    for accession_number in accession_numbers:
        fasta_string = get_fasta_from_ncbi(accession_number)
        if fasta_string:
            sequences = read_fasta_string(fasta_string)
            for header, sequence in sequences.items():
                gc_content = calculate_gc_content(sequence)
                length = len(sequence)
                statistic, p_value, observed_gc, expected_gc, sequence_length = perform_chisquare_test(sequence) # Добавил expected_gc

                print(f"Последовательность: {header} ({accession_number})")
                print(f"  Длина: {length}")
                print(f"  GC-состав: {gc_content:.2f}%")

                if statistic is not None and p_value is not None:
                    print(f"  Хи-квадрат: statistic = {statistic:.2f}, p-value = {p_value:.3f}")
                    if p_value < 0.05:
                        print("  GC-состав статистически значимо отличается от 50%.\n")
                    else:
                        print("  GC-состав статистически не отличается от 50%.\n")
                else:
                    print("  Недостаточно данных для проведения теста хи-квадрат.\n")

                visualize_gc_content(header, accession_number, observed_gc, expected_gc, statistic, p_value, sequence_length) # Передал expected_gc

        else:
            print(f"Не удалось получить последовательность для {accession_number}\n")

if __name__ == "__main__":
    main()