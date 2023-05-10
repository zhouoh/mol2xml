#include <iostream>
#include "tinyxml2.h"
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include "readmol.h"

using namespace std;
using namespace tinyxml2;

void traversingXML(XMLNode* node) {
    if(node == nullptr)
        return;
    
    if(node->ToDeclaration()) {
        auto declaration = dynamic_cast<XMLDeclaration*>(node);
        cout << "XML 声明，value=" << declaration->Value() << endl;
    }
    if(node->ToElement()) {
        auto element = dynamic_cast<XMLElement*>(node);
        cout << "XML 元素，name=" << element->Name() << ", value=" << element->Value() << endl;
        const XMLAttribute* attribute = element->FirstAttribute();
        while (attribute != nullptr) {
            cout << "\t属性 " << attribute->Name() << "=" << attribute->Value() << endl;
            attribute = attribute->Next();
        }
    }
    if(node->ToText()) {
        auto text = dynamic_cast<XMLText*>(node);
        cout << "XML 文本：" << text->Value() << endl;
    }
    if(node->ToComment()) {
        auto comment = dynamic_cast<XMLComment*>(node);
        cout << "XML 注释：" << comment->Value() << endl;
    }
    if(node->ToUnknown()) {
        auto unknown = dynamic_cast<XMLUnknown*>(node);
        cout << "XML 未知：" << unknown->Value() << endl;
    }
    if(node->ToDocument()) {
        auto document = dynamic_cast<XMLDocument*>(node);
        cout << "XML 文档：" << document->ErrorName() << endl;
    }
    
    if(node->NoChildren()) {
        return;
    }
    
    XMLNode* child = node->FirstChild();
    while(child != nullptr) {
        traversingXML(child);
        child = child->NextSibling();
    }
}

int main(int argc, const char * argv[]) {
    
    XMLDocument xmlDocument;
    XMLError error = xmlDocument.LoadFile("/Users/zhouoh/Downloads/scenes/cbox.xml");
    if(error != XML_SUCCESS) {
        std::cout << "读取 xml 失败：" << xmlDocument.ErrorStr() << endl;
        return EXIT_FAILURE;
    }
    
    traversingXML(&xmlDocument);
    return EXIT_SUCCESS;
}


