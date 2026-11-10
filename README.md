# Scientific Apps

A static GitHub Pages website containing a homepage and the latest FWI Radiation Lab. All calculations run in the browser. No paid server, API key or build step is required. The GitHub Free route uses a public repository and a public website.

## Publish on GitHub Pages

1. Create a **public** repository named `scientific-apps` in your own GitHub account. If you prefer your main user website instead, name it `<your-username>.github.io`.
2. Extract the ZIP and upload the contents of this folder into the repository root. `index.html`, `apps.json` and the `apps` folder must appear directly at the root, rather than inside another enclosing folder.
3. Commit the files to the `main` branch.
4. Open **Settings → Pages**. Under **Build and deployment**, choose **Deploy from a branch**, select **main** and **/(root)**, then click **Save**.
5. Wait for GitHub's Pages deployment to complete. The **Visit site** link in Settings → Pages is the authoritative published URL. GitHub notes that publishing can take up to ten minutes.

A project repository named `scientific-apps` normally uses `https://<your-username>.github.io/scientific-apps/`. A user repository named `<your-username>.github.io` uses `https://<your-username>.github.io/`. These are format examples, not links to an already published site.

Official instructions: https://docs.github.com/en/pages/getting-started-with-github-pages/creating-a-github-pages-site

## Add another app

1. Put its HTML/CSS/JavaScript under `apps/<app-slug>/`, including an `index.html`.
2. Add a record to `apps.json` with `id`, `name`, `category`, `description`, `path` and `tags`. Use a relative path such as `apps/my-tool/`.
3. Commit both the app and the catalog change. The homepage creates its card and category filter from `apps.json`.

Use relative asset paths so apps work beneath a GitHub project URL. For an app that requires Python or another backend, this website can link to an externally hosted instance; the current catalog renderer accepts locally hosted `apps/<slug>/` paths. Adapt its path validation deliberately before adding external links. The FWI app has a radiation-pattern preview. Other catalog entries receive a generic title-initials preview until an app-specific preview is added.

## Preview locally

Run this from the repository root:

    python3 -m http.server 8000

Then open http://localhost:8000. Use a server rather than opening index.html directly because the homepage reads apps.json with fetch.

## Verify the scientific app

    node checks/physics.cjs
    node checks/parameterizations.cjs
    node checks/review_regressions.cjs
    python3 checks/born_independent_check.py apps/fwi-radiation/physics.js

Python's independent check requires numpy. These dependencies are needed only for testing, not for hosting or using the app. SCIENTIFIC_CHANGES.md records the scientific corrections and remaining limits. SOURCE_REVISION.txt identifies the app revision.

## Files to publish

The archive contains static source, scientific checks and documentation. It excludes ChatGPT hosting settings, credentials, Git history, uploaded user images and private reports. No license is selected automatically; add an appropriate license before treating the source as licensed for reuse.
