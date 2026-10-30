# geodesy's website

The page for the geodesy library, live at https://geodesy-package.syncster.dev.

Author: Alexius Academia · Started 2026-10-31

A static page: no build step, no dependencies. Netlify (team syncster.dev)
deploys this folder by itself whenever `main` is pushed to GitHub.

```
index.html             # the page
styles.css             # the other libraries' design, in crimson
assets/route.png, scale.png   # the pictures examples/ draws with chart (copied from ../docs)
assets/share.png       # the picture shown when a link to the site is shared (1200 x 630)
assets/favicon.svg
netlify.toml           # no build; HTML isn't cached
```

To look at it locally: `python3 -m http.server 8080 -d website`, then
http://localhost:8080/.

The pictures come from `tessel run examples` in the library folder, which
writes them to docs/; copy them here. share.png is route.png cropped to
1372 x 720 and scaled to 1200 x 630.
