import os
import requests
import google.generativeai as genai

genai.configure(api_key=os.environ["GEMINI_API_KEY"])
model = genai.GenerativeModel('gemini-1.5-flash')

repo = os.environ["REPO_NAME"]
pr_number = os.environ["PR_NUMBER"]
token = os.environ["GITHUB_TOKEN"]

headers = {
    "Accept": "application/vnd.github.v3.diff",
    "Authorization": f"Bearer {token}"
}
url = f"https://api.github.com/repos/{repo}/pulls/{pr_number}"
response = requests.get(url, headers=headers)
diff_text = response.text

if not diff_text.strip():
    print("Diff пуст. Пропускаем ревью.")
    exit(0)

# Промпт адаптирован под проект AGO (оптимизация, шаблоны, OpenMP)
prompt = f"""
Ты Senior C++ Developer и эксперт по вычислительной математике. Проведи строгое code review для следующего git diff.
Проект (AGO) занимается решением задач глобальной оптимизации, активно использует C++ шаблоны (templates), многомерные массивы и многопоточность (OpenMP).

При ревью обрати особое внимание на:
1. Вычислительную эффективность (отсутствие лишних аллокаций в циклах, правильная передача по константной ссылке).
2. Корректность работы с шаблонами классов (например, TPoint, TInterval).
3. Безопасность многопоточности (условия гонки при использовании OpenMP).
4. Правильность работы с памятью и std::vector/std::array.
5. Соответствие стандартам современного C++ (C++17).

Выведи только конкретные замечания и куски кода с возможными исправлениями.

Diff:
```diff
{diff_text}
"""