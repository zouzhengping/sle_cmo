# GitHub 上传指南

## 方法一：使用命令行（推荐）

### 步骤 1: 在 GitHub 上创建新仓库

1. 登录 GitHub (https://github.com)
2. 点击右上角的 "+" 号，选择 "New repository"
3. 填写仓库信息：
   - **Repository name**: 例如 `sle_tscm_signature_analysis`
   - **Description**: "Analysis code for TSCM signature identification in SLE"
   - **Visibility**: 选择 Public 或 Private
   - **不要**勾选 "Initialize this repository with a README"（因为我们已经有了）
4. 点击 "Create repository"

### 步骤 2: 初始化本地 Git 仓库

在终端中执行以下命令：

```bash
# 进入代码文件夹
cd /home/ug2217/projects/sle_tscm_bulk/manuscript_code_for_github

# 初始化 Git 仓库
git init

# 添加所有文件
git add .

# 创建第一次提交
git commit -m "Initial commit: TSCM signature analysis code for SLE manuscript"
```

### 步骤 3: 连接到 GitHub 并推送

```bash
# 添加远程仓库（将 YOUR_USERNAME 和 REPO_NAME 替换为你的实际信息）
git remote add origin https://github.com/YOUR_USERNAME/REPO_NAME.git

# 或者使用 SSH（如果你配置了 SSH key）
# git remote add origin git@github.com:YOUR_USERNAME/REPO_NAME.git

# 推送代码到 GitHub
git branch -M main
git push -u origin main
```

**注意**: 如果这是第一次推送，GitHub 可能会要求你输入用户名和密码（或 Personal Access Token）。

---

## 方法二：使用 GitHub Desktop（图形界面）

如果你更喜欢使用图形界面：

1. 下载并安装 GitHub Desktop: https://desktop.github.com/
2. 打开 GitHub Desktop
3. 点击 "File" -> "Add Local Repository"
4. 选择文件夹: `/home/ug2217/projects/sle_tscm_bulk/manuscript_code_for_github`
5. 点击 "Publish repository" 按钮
6. 填写仓库名称和描述
7. 选择 Public 或 Private
8. 点击 "Publish Repository"

---

## 方法三：使用 GitHub Web 界面上传

1. 在 GitHub 上创建新仓库（步骤同方法一的步骤1）
2. 创建仓库后，GitHub 会显示上传文件的选项
3. 点击 "uploading an existing file"
4. 将 `manuscript_code_for_github` 文件夹中的所有文件拖拽到上传区域
5. 在页面底部填写提交信息，点击 "Commit changes"

---

## 常见问题

### 问题 1: 需要身份验证

如果推送时要求输入密码，你需要使用 **Personal Access Token** 而不是密码：

1. 登录 GitHub
2. 点击右上角头像 -> Settings
3. 左侧菜单选择 "Developer settings"
4. 选择 "Personal access tokens" -> "Tokens (classic)"
5. 点击 "Generate new token (classic)"
6. 设置权限（至少勾选 `repo`）
7. 生成后复制 token
8. 推送时，用户名输入你的 GitHub 用户名，密码输入刚才复制的 token

### 问题 2: 文件太大

如果某些文件太大无法上传，检查 `.gitignore` 文件是否已正确配置排除大文件。

### 问题 3: 更新代码

如果之后需要更新代码：

```bash
cd /home/ug2217/projects/sle_tscm_bulk/manuscript_code_for_github
git add .
git commit -m "Update: 描述你的更改"
git push
```

---

## 快速命令（复制粘贴）

将以下命令中的 `YOUR_USERNAME` 和 `REPO_NAME` 替换后，可以直接复制执行：

```bash
cd /home/ug2217/projects/sle_tscm_bulk/manuscript_code_for_github
git init
git add .
git commit -m "Initial commit: TSCM signature analysis code for SLE manuscript"
git branch -M main
git remote add origin https://github.com/YOUR_USERNAME/REPO_NAME.git
git push -u origin main
```

