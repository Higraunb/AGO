import os
import sys
import requests
from google import genai 

# Инициализация нового клиента
client = genai.Client(api_key=os.environ["GEMINI_API_KEY"])

repo = os.environ["REPO_NAME"]
pr_number = os.environ["PR_NUMBER"]
token = os.environ["GITHUB_TOKEN"]

print(f"Запуск ревью для PR #{pr_number} в репозитории {repo}...", flush=True)

# 1. Получаем diff
headers = {
    "Accept": "application/vnd.github.v3.diff",
    "Authorization": f"Bearer {token}",
    "X-GitHub-Api-Version": "2022-11-28"
}
url = f"https://api.github.com/repos/{repo}/pulls/{pr_number}"
response = requests.get(url, headers=headers)

if response.status_code != 200:
    print(f"Ошибка получения diff от GitHub: Код {response.status_code}\nОтвет: {response.text}", flush=True)
    sys.exit(1)

diff_text = response.text

if not diff_text.strip():
    print("Diff пуст (нет изменений в коде). Пропускаем ревью.", flush=True)
    sys.exit(0)

# 2. Промпт
prompt = f"""
Ты Senior C++ Developer и эксперт по вычислительной математике. Проведи строгое code review для следующего git diff.
Проект (AGO) занимается решением задач глобальной оптимизации, активно использует C++ шаблоны (templates), многомерные массивы и многопоточность (OpenMP).

При ревью обрати особое внимание на:
1. Вычислительную эффективность (отсутствие лишних аллокаций в циклах, правильная передача по константной ссылке).
2. Корректность работы с шаблонами классов.
3. Безопасность многопоточности (OpenMP).
4. Правильность работы с памятью и std::vector/std::array.
5. Соответствие стандартам современного C++ (C++17).

Выведи только конкретные замечания и куски кода с возможными исправлениями.

Diff:
```diff
{diff_text}
"""

print("Отправка кода на анализ в Gemini...", flush=True)

try:
    ai_response = client.models.generate_content(
        model='gemini-2.5-flash',
        contents=prompt
    )
    review_comment = ai_response.text
except Exception as e:
    print(f"Ошибка при вызове Gemini API: {e}", flush=True)
    sys.exit(1)

print("Публикация ответа на GitHub...", flush=True)

comment_url = f"https://api.github.com/repos/{repo}/issues/{pr_number}/comments"
comment_headers = {
"Accept": "application/vnd.github.v3+json",
"Authorization": f"Bearer {token}",
"X-GitHub-Api-Version": "2022-11-28"
}

resp = requests.post(comment_url, headers=comment_headers, json={"body": f"### 🤖 Gemini AI Code Review\n\n{review_comment}"})

print(f"Статус отправки комментария: {resp.status_code}", flush=True)
if resp.status_code != 201:
  print(f"Ошибка публикации от GitHub: {resp.text}", flush=True)
else:
  print("Успешно!", flush=True)