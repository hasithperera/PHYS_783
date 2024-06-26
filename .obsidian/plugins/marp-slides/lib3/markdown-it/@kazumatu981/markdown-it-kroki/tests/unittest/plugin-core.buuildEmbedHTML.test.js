const md = require('markdown-it');
const expect = require('chai').expect;
const { JSDOM } = require('jsdom');
const { MarkdownItKrokiCore } = require('../../lib/plugin-core');
const { encode } = require('../../lib/diagram-encoder');

describe('# [unit-test] plugin-core.js', () => {
    describe('## method: buuildEmbedHTML', () => {
        describe('### langAndAlt.language', () => {
            [null, undefined, ''].forEach((test) => {
                it(`* when langAndAlt.language is null or empty, throws error. testcase:${test}`, () => {
                    const diagramCode = '@startuml\nBob -> Alice : hello\n @enduml';
                    const plugin = new MarkdownItKrokiCore(new md()).setOptions();
                    const testFunc = () => {
                        plugin.use();
                        const _ = plugin.buildEmbedHTML(
                            { language: test, alt: '' }, diagramCode);
                    };
                    expect(testFunc).to.throw();
                });
            });
            it('* language embeded in to url', () => {
                const test = 'plantuml';

                const diagramCode = '@startuml\nBob -> Alice : hello\n @enduml';

                // build embed HTML
                const plugin = new MarkdownItKrokiCore(new md()).setOptions();
                plugin.use();
                const html = plugin.buildEmbedHTML({ language: test, alt: '' }, diagramCode);

                // parse dom
                const dom = new JSDOM(html);
                const imgTag = dom.window.document.getElementsByTagName("embed")[0];

                // get url attribute
                const url = imgTag.getAttribute('src');

                expect(/\/plantuml\//.test(url)).to.true;
            });
        });
        describe('### langAndAlt.alt', () => {
            [null, undefined, ''].forEach((test) => {
                it(`* when langAndAlt.alt is null or empty, no alt attribute. testcase:${test}`, () => {
                    const diagramCode = '@startuml\nBob -> Alice : hello\n @enduml';

                    // prepair
                    const plugin = new MarkdownItKrokiCore(new md()).setOptions();
                    plugin.use();

                    // render
                    const html = plugin.buildEmbedHTML(
                        { language: 'plantuml', alt: test }, diagramCode);
                    // parse dom
                    const dom = new JSDOM(html);
                    const imgTag = dom.window.document.getElementsByTagName("embed")[0];

                    expect(imgTag.hasAttribute('alt')).to.false;
                });
            });
            it('* embeded altText', () => {
                const expected = "this is test Text";
                const diagramCode = '@startuml\nBob -> Alice : hello\n @enduml';

                // prepair
                const plugin = new MarkdownItKrokiCore(new md()).setOptions();
                plugin.use();

                // render
                const html = plugin.buildEmbedHTML(
                    { language: 'plantuml', alt: expected }, diagramCode);
                // parse dom
                const dom = new JSDOM(html);
                const imgTag = dom.window.document.getElementsByTagName("embed")[0];

                expect(imgTag.getAttribute('title')).to.equal(expected);

            });
        });
        describe('### diagramCode', () => {
            [null, undefined, ''].forEach((test) => {
                it(`* when diagramCode is null or empty, throws error. testcase:${test}`, () => {
                    const plugin = new MarkdownItKrokiCore(new md()).setOptions();
                    const testFunc = () => {
                        plugin.use();
                        const _ = plugin.buildEmbedHTML(
                            { language: 'plantuml', alt: '' }, test);
                    };
                    expect(testFunc).to.throw();
                });
            });
            it('* encoded diagram must be embed to url on <img src=\'....\' ', () => {
                const test = '@startuml\nBob -> Alice : hello\n @enduml';

                const expected = encode(test);

                // build embed HTML
                const plugin = new MarkdownItKrokiCore(new md()).setOptions();
                plugin.use();
                const html = plugin.buildEmbedHTML({ language: 'plantuml', alt: '' }, test);

                // parse dom
                const dom = new JSDOM(html);
                const imgTag = dom.window.document.getElementsByTagName("embed")[0];

                // get url attribute
                const url = imgTag.getAttribute('src');

                expect(url.endsWith(expected)).to.true;
            });
            it('* <img> is surounded by <marp-auto-scaling> on to be used form marp-it', () => {
                const test = '@startuml\nBob -> Alice : hello\n @enduml';

                const markdownIt = new md()
                markdownIt['marpit'] = { someObject: 'is implemented' };
                // build embed HTML
                const plugin = new MarkdownItKrokiCore(markdownIt).setOptions();
                plugin.use();
                const html = plugin.buildEmbedHTML({ language: 'plantuml', alt: '' }, test);

                // parse dom
                const dom = new JSDOM(html);
                const imgTag = dom.window.document.getElementsByTagName("embed")[0];
                const marpAutoScaling = dom.window.document.getElementsByTagName("marp-auto-scaling")[0];

                expect(imgTag.isSameNode(marpAutoScaling.firstChild)).to.be.true;
            });
            it('* <img> is surounded by <marp-auto-scaling> on not to be used form marp-it', () => {
                const test = '@startuml\nBob -> Alice : hello\n @enduml';

                const markdownIt = new md()
                // build embed HTML
                const plugin = new MarkdownItKrokiCore(markdownIt).setOptions();
                plugin.use();
                const html = plugin.buildEmbedHTML({ language: 'plantuml', alt: '' }, test);

                // parse dom
                const dom = new JSDOM(html);
                const marpAutoScaling = dom.window.document.getElementsByTagName("marp-auto-scaling");

                expect(marpAutoScaling.length).to.equal(0);
            })

        });
    });
});